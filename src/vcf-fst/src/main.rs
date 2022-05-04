use fst::set::{Stream, StreamBuilder};
use std::path::{PathBuf, Path};
use std::error::Error;

use std::fs::File;
use fst::{IntoStreamer, SetBuilder, Set};
//use fst::Streamer;
use fst::automaton::{Automaton, Str, StartsWith};
//use fst::automaton::Subsequence;
//use memmap::Mmap;

use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use rayon; 
use pedigree_sims::io::vcf::{
    self,
    SampleTag,
    reader::{VCFReader, VCFAsyncReader, VCFPanelReader},
};

pub struct VCFIndexer<'a>{
    reader              : VCFReader<'a>,                      // VCFReader.
    builder             : SetBuilder<std::io::BufWriter<File>>,   // FST SetBuilder
    frq_builder        : SetBuilder<std::io::BufWriter<File>>,
    //samples             : Vec<io::SampleTag>,                   // Sorter vector of sampletags, sorted across Id.
    samples             : BTreeMap<SampleTag, String>,
    counter             : usize,                                  // Line counter.
    coordinate_buffer   : Vec<u8>,                                // Coordinate of the current line
    previous_coordinate : Vec<u8>,                                // Coordinate of the previous line
    genotypes_buffer    : Vec<u8>,                                // Genotypes  of the current line.
    genotype_keys       : Vec<Vec<u8>>,                           // Values that should be inserted into the set builder.
    frequency_keys      : Vec<String>
}      

impl<'a> VCFIndexer<'a> {
    pub  fn new(vcf: &Path, output_file: &str, samples: BTreeMap<SampleTag, String>, threads: usize) -> Result<Self, Box<dyn Error>> {
        // ----------------------- Initialize Writer
        let gen_writer  = std::io::BufWriter::new(File::create(format!("{output_file}.fst"))?);
        let builder =  SetBuilder::new(gen_writer)?;

        let frq_writer  = std::io::BufWriter::new(File::create(format!("{output_file}.frq.fst"))?);
        let frq_builder =  SetBuilder::new(frq_writer)?;

        // ----------------------- Initialize Reader.
        let reader  = VCFReader::new(vcf, threads)?;
        // ----------------------- Preprocess samples ID and Index.
        //let samples = Self::parse_samples(&reader);
        Ok(Self {
            reader,
            builder,
            frq_builder,
            samples,
            counter            : 0,
            coordinate_buffer  : Vec::new(),
            previous_coordinate: Vec::new(),
            genotypes_buffer   : Vec::new(),
            genotype_keys      : Vec::new(),
            frequency_keys     : Vec::new(),
        })
    }

    pub  unsafe fn build_fst(&mut self) -> Result<(), Box<dyn Error>> {
        'line: while self.reader.has_data_left()? {

            self.fill_current_position_buffer()?; // Parse the current_position
            if self.insert_previous_keys()? {     // Check if the previous line had the same coordinates. Do NOT insert keys of the previous line if it is the case.
                continue 'line
            }
            self.save_previous_line();            // Save the current coordinate for the next iteration
    

            // Go to INFO field.
            self.skip_fields(5)?;                                      // 6 
            let info = self.reader.next_field()?.split(';').collect::<Vec<&str>>();

            // Check if Bi-Allelic and skip line if not.
            if info.iter().any(|&field| field == "MULTI_ALLELIC") {
                //println!("Not Biallelic. Skipping.");
                self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
                self.coordinate_buffer.clear();
                continue 'line
            }

            // Get the mutation type from INFO...
            let vtype = info.iter()
                .find(|&&field| field.starts_with("VT=")).unwrap()
                .split('=')
                .collect::<Vec<&str>>()[1];

            // ...and skip line if this is not an SNP.
            if vtype != "SNP" {
                //println!("Not an SNP. Skipping.");
                self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
                self.coordinate_buffer.clear();
                continue 'line
            }
            // Extract population allele frequency.
            let mut pop_afs = info.iter()
                .filter(|&&field| field.contains("_AF"))
                .map(|field| field.split("_AF=").collect())
                .collect::<Vec<Vec<&str>>>();          
            pop_afs.sort();

            for pop_af in pop_afs.iter() {
                let pop_af_key = format!("{}{} {}", std::str::from_utf8(&self.coordinate_buffer)?, pop_af[0], pop_af[1]);
                self.frequency_keys.push(pop_af_key);
            }


            self.skip_fields(1)?;                 // Skip INFO field
            self.fill_genotypes_buffer()?;      
            self.print_progress(50000)?;              
            self.parse_keys()?;
            //println!("{:?}", self.keys);
            self.genotypes_buffer.clear();
            self.coordinate_buffer.clear();
        }

        self.insert_previous_keys()?;
        Ok(())
    }

    pub fn finish_build(self) -> Result<(), Box<dyn Error>> {
        self.builder.finish()?;
        self.frq_builder.finish()?;
        Ok(())
    }

    fn parse_samples(reader: &VCFAsyncReader) -> Vec<SampleTag> {
        // Extract and sort samples 
        let samples = reader.samples();
        let samples = &samples[9..];
        // ----------------------- Convert to struct to keep id.
        let mut struct_samples = Vec::new();
        for (i, sample) in samples.iter().enumerate(){
            let sample = SampleTag::new(sample, Some(i));
            struct_samples.push(sample);
        }
        // ----------------------- Sort samples by ID.
        struct_samples.sort();
        struct_samples
    }

     fn fill_current_position_buffer(&mut self) -> Result<(), Box<dyn Error>> {
        // ---------------------------- Get current chromosome and position.
        let chr_bytes = self.reader.source.read_until(b'\t', &mut self.coordinate_buffer)?; // Add chromosome
        self.add_field_separator(b' ')?;
        let pos_bytes = self.reader.source.read_until(b'\t', &mut self.coordinate_buffer)?; // Add position.
        self.add_field_separator(b' ')?;
        for _ in 0..(10-pos_bytes) {                 // Add 9 leading zeroes for padding (this ensure our key is sorted.)
            self.coordinate_buffer.insert(chr_bytes, b'0')   
        }
        Ok(())
    }

     fn add_field_separator(&mut self, field: u8) -> Result<(), Box<dyn Error>> {
        self.coordinate_buffer.pop();
        self.coordinate_buffer.push(field);
        Ok(())
    }
    
    fn duplicate_line(&self) -> bool {
        self.coordinate_buffer == self.previous_coordinate
    }


     fn insert_previous_keys(&mut self) -> Result<bool, Box<dyn Error>> {
        if self.duplicate_line() {
            //println!("Duplicate line at : [{}]. Skipping.", std::str::from_utf8(&self.previous_coordinate)?);
            self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
            self.coordinate_buffer.clear();
            Ok(true)
        } else {
            //println!("Inserting keys.");
            self.insert_keys()?;       // If the line is different, we can then insert the keys of the previous line.
            Ok(false)
        }
    }

     fn clear_buffers(&mut self) -> Result<(), Box<dyn Error>> {
        self.reader.source.read_until(b'\n', &mut Vec::new())?;
        self.genotype_keys.clear();
        self.frequency_keys.clear();
        Ok(())
    }

     fn insert_keys(&mut self) -> Result<(), Box<dyn Error>> {
        for gen in self.genotype_keys.iter() {
            self.builder.insert(gen)?;
        }
        self.genotype_keys.clear();

        for frq in self.frequency_keys.iter() {
            self.frq_builder.insert(frq)?;
        } 
        self.frequency_keys.clear();
        Ok(())
    }

    fn save_previous_line(&mut self) {
        self.previous_coordinate = self.coordinate_buffer.clone();
    }

     fn skip_fields(&mut self, n: usize) -> Result<(), Box<dyn Error>> {
        for _ in 0..n {
            self.reader.source.read_until(b'\t', &mut Vec::new())?;
        }    
        Ok(())
    }

     fn fill_genotypes_buffer(&mut self) -> Result<(), Box<dyn Error>> {   
        self.reader.source.read_until(b'\n', &mut self.genotypes_buffer)?;
        Ok(())
    }

    fn print_progress(&mut self, every_n: usize) -> Result<(), Box<dyn Error>> {
        if self.counter % every_n == 0 {
            println!("{: >9} {:?}", self.counter, std::str::from_utf8(&self.coordinate_buffer)?);
        }
        self.counter+= 1;
        Ok(())
    }

    unsafe fn parse_keys(&mut self) -> Result<(), Box<dyn Error>> {
        for (sample_tag, sample_pop) in self.samples.iter() {            

            // Complete value to insert.
            let mut key = self.coordinate_buffer.clone();
            key.append(sample_tag.id().clone().as_mut_vec());  // UNSAFE: add sampleID
            key.push(b' ');

            // Extract genotype info for thesample and add it to the key.
            let geno_idx = sample_tag.idx().unwrap()*4;                // x4, since each field+separator is four characters long. (e.g.: '0|0 ')
            key.push(self.genotypes_buffer[geno_idx  ]); // haplo1
            key.push(self.genotypes_buffer[geno_idx+2]); // haplo2

            //println!("{}", std::str::from_utf8(&key)?);
            self.genotype_keys.push(key);
        }
        Ok(())
    }
}




fn main() -> Result<(), Box<dyn Error>> {
    // ---------- Command line arguments
    let cli_threads = 22;
    let decompression_threads= 0;
    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b/");
    let panel = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b/integrated_call_samples_v3.20130502.ALL.panel");
    let user_defined_subset: Option<Vec<&str>>  = Some(vec!["EUR", "AFR"]);
    // ------------------------------------

    println!("Fetching input VCF files in {}", &data_dir.to_str().unwrap());
    let mut input_vcf_paths = vcf::get_input_vcfs(&data_dir).unwrap();
    input_vcf_paths.sort();


    let mut panel = VCFPanelReader::new(panel.as_path())?;
    panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;

    // User defined subset-arguments.
    match &user_defined_subset {
        Some(subset) => panel.subset_panel(&subset.clone()),
        None         => (),
    }

    // --------------------- Set ThreadPool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cli_threads)
        .build()
        .unwrap();


    pool.scope(|scope|{
        for vcf in input_vcf_paths.iter(){
            scope.spawn(|_| {
                // -------------------------- Format output filename
                let file_stem = vcf.as_path()
                    .file_name()
                    .unwrap()
                    .to_str().unwrap()
                    .replace(".vcf.gz", "")
                    .replace(".vcf", "");

                let pop_tag = match &user_defined_subset {
                    Some(subset) => format!("-{}", subset.join("-")),
                    None => "".to_string(),
                };
                let output_path =format!("tests/test-data/fst/{}{}", file_stem, pop_tag);
                println!("{output_path:?}");
                let mut setbuilder = VCFIndexer::new(vcf, &output_path, panel.into_transposed_btreemap(), decompression_threads).unwrap();
                unsafe {setbuilder.build_fst().unwrap();}
                setbuilder.finish_build().unwrap();
            });
        }
    });




    // Bincode strategy

    //serialize(&input_vcf_paths, &tpool).unwrap();

    // FST Strategy
    println!("Generating FST index...");
    //sort_vcf_map(&input_vcf_paths, &tpool).unwrap();
    //vcf_fst_index(&input_vcf_paths[0], 4)?;
    println!("Done!");

    //// Memory Map strategy.
    //println!("Opening memory map");
    //let mmap = unsafe { Mmap::map(&File::open("tests/fst/g1k-map.fst").unwrap()).unwrap() };
    //let map = Map::new(mmap).unwrap();

    //// In RAM Strategy
    println!("Reading in memory.");
    let mut file_handle = File::open("tests/test-data/fst/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-EUR-AFR.fst").unwrap();
    let mut bytes = vec![];
    std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
    //let map = Map::new(bytes).unwrap();

    let set = Set::new(bytes).unwrap();

    let samples = ["HG01067", "NA20885", "HG00267", "HG00236", "NA20356", "NA19346"];
    let samples =  ["HG00096", "HG00264", "HG01618", "HG01920"]; //EUR + AFR

    for sample in samples.iter() {
        println!("Searching {sample}");
        //000061543
        //000030524
        let regex = format!("1 000055164 {sample}");
        //let matcher = Subsequence::new(&regex);
        let matcher = Str::new(&regex)
            .starts_with();
            //.intersection(Subsequence::new(&sample));
        let stream = set.search(&matcher).into_stream();

        println!("  - Populating vector");
        //let mut kvs = vec![];
        let kvs = stream.into_strs().unwrap();
        //let stream = map.search(&matcher).into_stream();
        //let kvs = stream.into_str_vec().unwrap();
        
        println!("{kvs:?}");
        //while let Some((k, v)) = stream.next() {
        //    println!("{:?} {}", k, v);
        //    let record = (String::from_utf8(k.to_vec()).unwrap(), v);
        //    println!("{:?}", &record);
        //    kvs.push(record);
        //    break
        //}
        //println!("{kvs:?}\n");
    }

    Ok(())
}
