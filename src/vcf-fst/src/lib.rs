use std::{
    path::Path,
    error::Error,
    fs::File,
    io::BufRead,
    collections::BTreeMap,
};

use fst::{
    SetBuilder,
};

use parser::{self};
use pedigree_sims::io::vcf::{
    self,
    SampleTag,
    reader::{VCFReader, VCFPanelReader},
};

use log::{info};
use rayon::{self}; 

pub struct VCFIndexer<'a>{
    reader              : VCFReader<'a>,                          // VCFReader.
    builder             : SetBuilder<std::io::BufWriter<File>>,   // FST SetBuilder for genotypes
    frq_builder         : SetBuilder<std::io::BufWriter<File>>,   // FST SetBuilder for allele_frequencies.
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

        let frq_writer  = std::io::BufWriter::new(File::create(format!("{output_file}.fst.frq"))?);
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
        for sample_tag in self.samples.keys() {            

            // Complete value to insert.
            let mut key = self.coordinate_buffer.clone();
            key.append(sample_tag.id().clone().as_mut_vec());  // UNSAFE: add sampleID
            key.push(b' ');

            // Extract genotype info for thesample and add it to the key.
            let geno_idx = sample_tag.idx().unwrap()*4;                // x4, since each field+separator is four characters long. (e.g.: '0|0 ')
            key.push(self.genotypes_buffer[geno_idx  ]); // haplo1
            key.push(self.genotypes_buffer[geno_idx+2]); // haplo2

            self.genotype_keys.push(key);
        }
        Ok(())
    }
}


pub fn run(
    fst_cli: &parser::VCFFst,
) -> Result<(), Box<dyn Error>> {

    info!("Fetching input VCF files in {}", &fst_cli.vcf_dir.to_str().unwrap());
    let mut input_vcf_paths = vcf::get_input_vcfs(&fst_cli.vcf_dir).unwrap();
    input_vcf_paths.sort();

    // --------------------- Fetch the input panel.
    let input_panel_path = match fst_cli.panel.clone() {
        Some(path) => path,
        None => vcf::fetch_input_panel(&fst_cli.vcf_dir)?,
    };

    let mut panel = VCFPanelReader::new(&input_panel_path)?;
    panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;

    // User defined subset-arguments.
    match &fst_cli.pop_subset {
        Some(subset) => panel.subset_panel(subset),
        None         => (),
    }

    // --------------------- Set ThreadPool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(fst_cli.threads)
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

                let pop_tag = match &fst_cli.pop_subset {
                    Some(subset) => format!("-{}", subset.join("-")),
                    None => "".to_string(),
                };
                let output_path =format!("{}/{}{}", fst_cli.output_dir.to_str().unwrap(), file_stem, pop_tag);
                println!("{output_path:?}");
                let mut setbuilder = VCFIndexer::new(vcf, &output_path, panel.into_transposed_btreemap(), fst_cli.decompression_threads).unwrap();
                unsafe {setbuilder.build_fst().unwrap();}
                setbuilder.finish_build().unwrap();
            });
        }
    });


    let output_panel_path = format!("{}/{}",
        fst_cli.output_dir.to_str().unwrap(),
        input_panel_path.file_name().unwrap().to_str().unwrap());
    panel.copy_from_source(Path::new(&output_panel_path)).unwrap();
    Ok(())
}
