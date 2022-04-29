use pedigree_sims::io;
use std::path::{PathBuf, Path};
use std::error::Error;

use std::fs::File;
use fst::{IntoStreamer, SetBuilder, Set};
//use fst::Streamer;
use fst::automaton::{Automaton, Str};
//use fst::automaton::Subsequence;
//use memmap::Mmap;

use std::collections::BTreeMap;

use tokio::io::AsyncBufReadExt;
use tokio;

pub struct VCFIndexer{
    reader              : io::VCFAsyncReader,                   // VCFReader.
    builder             : SetBuilder<std::io::BufWriter<File>>, // FST SetBuilder
    //samples             : Vec<io::SampleTag>,                   // Sorter vector of sampletags, sorted across Id.
    samples             : BTreeMap<io::SampleTag, String>,
    counter             : usize,                                // Line counter.
    coordinate_buffer   : Vec<u8>,                              // Coordinate of the current line
    previous_coordinate : Vec<u8>,                              // Coordinate of the previous line
    genotypes_buffer    : Vec<u8>,                              // Genotypes  of the current line.
    keys                : Vec<Vec<u8>>,                         // Values that should be inserted into the set builder.
}      

impl VCFIndexer {
    pub async fn new(vcf: &Path, output_file: &str, samples: BTreeMap<io::SampleTag, String>, threads: usize) -> Result<Self, Box<dyn Error>> {
        // ----------------------- Initialize Writer
        let writer  = std::io::BufWriter::new(File::create(output_file)?);
        let builder =  SetBuilder::new(writer)?;
        // ----------------------- Initialize Reader.
        let reader  = io::VCFAsyncReader::new(vcf, threads).await?;
        // ----------------------- Preprocess samples ID and Index.
        //let samples = Self::parse_samples(&reader);
        Ok(Self {
            reader,
            builder,
            samples,
            counter            : 0,
            coordinate_buffer  : Vec::new(),
            previous_coordinate: Vec::new(),
            genotypes_buffer   : Vec::new(),
            keys               : Vec::new(),
        })
    }

    pub async unsafe fn build_fst(&mut self) -> Result<(), Box<dyn Error>> {
        'line: while self.reader.has_data_left().await? {

            self.fill_current_position_buffer().await?; // Parse the current_position
            if self.insert_previous_keys().await? {     // Check if the previous line had the same coordinates. Do NOT insert keys of the previous line if it is the case.
                continue 'line
            }
            self.save_previous_line();                  // Save the current coordinate for the next iteration
    
            self.skip_fields(7).await?;                 // Skip INFO field
            self.fill_genotypes_buffer().await?;      
            self.print_progress(50000)?;              
            self.parse_keys()?;
            //println!("{:?}", self.keys);
            self.genotypes_buffer.clear();
            self.coordinate_buffer.clear();
        }

        self.insert_previous_keys().await?;
        Ok(())
    }

    pub fn finish_build(self) -> Result<(), Box<dyn Error>> {
        self.builder.finish()?;
        Ok(())
    }

    fn parse_samples(reader: &io::VCFAsyncReader) -> Vec<io::SampleTag> {
        // Extract and sort samples 
        let samples = reader.samples();
        let samples = &samples[9..];
        // ----------------------- Convert to struct to keep id.
        let mut struct_samples = Vec::new();
        for (i, sample) in samples.iter().enumerate(){
            let sample = io::SampleTag::new(sample, i);
            struct_samples.push(sample);
        }
        // ----------------------- Sort samples by ID.
        struct_samples.sort();
        struct_samples
    }

    async fn fill_current_position_buffer(&mut self) -> Result<(), Box<dyn Error>> {
        // ---------------------------- Get current chromosome and position.
        let chr_bytes = self.reader.source.read_until(b'\t', &mut self.coordinate_buffer).await?; // Add chromosome
        self.add_field_separator(b' ').await?;
        let pos_bytes = self.reader.source.read_until(b'\t', &mut self.coordinate_buffer).await?; // Add position.
        self.add_field_separator(b' ').await?;
        for _ in 0..(10-pos_bytes) {                 // Add 9 leading zeroes for padding (this ensure our key is sorted.)
            self.coordinate_buffer.insert(chr_bytes, b'0')   
        }
        Ok(())
    }

    async fn add_field_separator(&mut self, field: u8) -> Result<(), Box<dyn Error>> {
        self.coordinate_buffer.pop();
        self.coordinate_buffer.push(field);
        Ok(())
    }
    
    fn duplicate_line(&self) -> bool {
        self.coordinate_buffer == self.previous_coordinate
    }


    async fn insert_previous_keys(&mut self) -> Result<bool, Box<dyn Error>> {
        if self.duplicate_line() {
            println!("Duplicate line at : [{}]. Skipping.", std::str::from_utf8(&self.previous_coordinate)?);
            self.clear_buffers().await?;     // If so, skip this lines, and don't insert the keys of the previous line.
            self.coordinate_buffer.clear();
            Ok(true)
        } else {
            self.insert_keys().await?;       // If the line is different, we can then insert the keys of the previous line.
            Ok(false)
        }
    }

    async fn clear_buffers(&mut self) -> Result<(), Box<dyn Error>> {
        self.reader.source.read_until(b'\n', &mut Vec::new()).await?;
        self.keys.clear();
        Ok(())
    }

    async fn insert_keys(&mut self) -> Result<(), Box<dyn Error>> {
        for key in self.keys.iter() {
            self.builder.insert(key)?;
        }
        self.keys.clear();
        Ok(())
    }

    fn save_previous_line(&mut self) {
        self.previous_coordinate = self.coordinate_buffer.clone();
    }

    async fn skip_fields(&mut self, n: usize) -> Result<(), Box<dyn Error>> {
        for _ in 0..n {
            self.reader.source.read_until(b'\t', &mut Vec::new()).await?;
        }    
        Ok(())
    }

    async fn fill_genotypes_buffer(&mut self) -> Result<(), Box<dyn Error>> {   
        self.reader.source.read_until(b'\n', &mut self.genotypes_buffer).await?;
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
            let geno_idx = sample_tag.idx()*4;                // x4, since each field+separator is four characters long. (e.g.: '0|0 ')
            key.push(self.genotypes_buffer[geno_idx  ]); // haplo1
            key.push(self.genotypes_buffer[geno_idx+2]); // haplo2

            //println!("{}", std::str::from_utf8(&key)?);
            self.keys.push(key);
        }
        Ok(())
    }
}



#[tokio::main(flavor = "multi_thread", worker_threads = 4)]
async fn main() -> Result<(), Box<dyn Error>> {

    
    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b-EUR-AFR-filtered");
    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b-first-50k-filtered/");
    let data_dir = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b-first-5000/");
    println!("Fetching input VCF files in {}", &data_dir.to_str().unwrap());
    let mut input_vcf_paths = io::get_input_vcfs(&data_dir).unwrap();
    input_vcf_paths.sort();


    let panel = PathBuf::from("tests/test-data/vcf/g1k-phase3-v5b/integrated_call_samples_v3.20130502.ALL.panel");
    let mut panel = io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path()).await?;

    // User defined subset-arguments.
    let user_defined_subset = Some(vec!["EUR", "AFR"]);
    match &user_defined_subset {
        Some(subset) => panel.subset_panel(&subset.clone()),
        None         => (),
    }


    for vcf in input_vcf_paths.iter(){
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
        let output_path =format!("tests/test-data/fst/{}{}.fst", file_stem, pop_tag);
        println!("{output_path:?}");
        let mut setbuilder = VCFIndexer::new(vcf, &output_path, panel.into_transposed_btreemap(), 4).await?;
        unsafe {setbuilder.build_fst().await?;}
        setbuilder.finish_build()?;
    }




    // Bincode strategy

    //serialize(&input_vcf_paths, &tpool).unwrap();

    // FST Strategy
    println!("Generating FST index...");
    //sort_vcf_map(&input_vcf_paths, &tpool).unwrap();
    //vcf_fst_index(&input_vcf_paths[0], 4).await?;
    println!("Done!");

    //// Memory Map strategy.
    //println!("Opening memory map");
    //let mmap = unsafe { Mmap::map(&File::open("tests/fst/g1k-map.fst").unwrap()).unwrap() };
    //let map = Map::new(mmap).unwrap();

    //// In RAM Strategy
    println!("Reading in memory.");
    let mut file_handle = File::open("tests/test-data/fst/g1k-set.fst").unwrap();
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
