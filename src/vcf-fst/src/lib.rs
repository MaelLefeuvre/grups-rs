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

/// Convert a `.vcf(.gz)` file into a Finite-State-Transducer Set.
/// See the following resources for a recap on the theory and applications of FST Sets :
/// - https://cs.nyu.edu/~mohri/pub/fla.pdf (DOI: 10.1007/978-3-540-39886-8_29) 
/// - https://blog.burntsushi.net/transducers/
/// 
/// # Behavior:
/// `VCFIndexer` will build to FST sets from a unique vcf file:
/// - The first set is identified by the `.fst` file extension, and contains the genotype information 
///   of each sample, at each coordinate.
///   - Fields (space-separated): <CHR>    <POS>    <SampleId>    <ALLELE1><ALLELE2>
/// 
/// - The second set is identified by the `.fst.frq` file extension and contains population allele frequencies
///   - Fields (space-separated): <CHR>    <POS>    <POP>    <FREQUENCY>
/// 
/// - <POS> fields are padded with up to 9 leading '0' characters, to ensure they are properly sorted.
/// 
/// `VCFIndexer` will automatically filter out unrelevant positions from the vcf file. This include:
///  - Duplicate coordinate entries
///  - Coordinates containing a `MULTI_ALLELIC` tag within the `INFO` field
///  - Indels (entries which are not of type `VT=SNP` are excluded from the set)
/// 
/// # Struct Fields
/// - `reader`             : inner `VCFReader`
/// - `builder`            : FST SetBuilder for genotypes (this is what constructs our `.fst` file)
/// - `frq_builder`        : FST SetBuilder for alleleq_frequencies (this is what constructs ou `.fst.frq` file)
/// - `samples`            : Transposed input panel definition map, sorted across `SampleTag.id()`
/// - `counter`            : Keeps track of the number of processed VCF lines.
/// - `coordinate_buffer`  : Coordinate of the current VCF line
/// - `previous_coordinate`: Coordinate of the previous VCF line
/// - `genotypes_buffer`   : Genotypes of the current VCF line
/// - `genotypes_keys`     : Temporary buffer of the values that should be inserted into the genotypes set builder
/// - `frequency_keys`     : Temporary buffer of the values that should be inserted into the frequency set builder
pub struct VCFIndexer<'a>{
    reader              : VCFReader<'a>,
    builder             : SetBuilder<std::io::BufWriter<File>>,
    frq_builder         : SetBuilder<std::io::BufWriter<File>>,
    samples             : BTreeMap<SampleTag, String>,
    counter             : usize,
    coordinate_buffer   : Vec<u8>,
    previous_coordinate : Vec<u8>,
    genotypes_buffer    : Vec<u8>,
    genotype_keys       : Vec<Vec<u8>>,
    frequency_keys      : Vec<String>
}      

impl<'a> VCFIndexer<'a> {
    /// Initialize a new VCFIndexer.
    /// # Arguments:
    /// - `vcf`    : path leading to the target `.vcf`|`.vcf.gz` file
    /// - `samples`: transposed input samples definition file (sorted across SampleTag.id())
    /// - `threads`: number of requested decompression threads (only relevant in the case of BGZF compressed `.vcf.gz` files)
    pub  fn new(vcf: &Path, output_file: &str, samples: BTreeMap<SampleTag, String>, threads: usize) -> Result<Self, Box<dyn Error>> {
        // ----------------------- Initialize Genotypes and Frequency Writers
        let gen_writer  = std::io::BufWriter::new(File::create(format!("{output_file}.fst"))?);
        let builder =  SetBuilder::new(gen_writer)?;

        let frq_writer  = std::io::BufWriter::new(File::create(format!("{output_file}.fst.frq"))?);
        let frq_builder =  SetBuilder::new(frq_writer)?;

        // ----------------------- Initialize VCFReader.
        let reader  = VCFReader::new(vcf, threads)?;
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

    /// Read through the contents of `self.reader` and dynamically build our `.fst` and `.fst.frq` sets.
    /// # Safety
    /// This function is unsafe because parsing keys requires constructing &mut Vec without utf8 validation.
    /// If this constraint is violated, using the original String after dropping the &mut Vec may violate memory safety
    /// as Rust assumes that Strings are valid UTF-8
    pub unsafe fn build_fst(&mut self) -> Result<(), Box<dyn Error>> {
        'line: while self.reader.has_data_left()? {
            // ---- Parse the current position 
            self.fill_current_position_buffer()?;

            // ---- Flush the previous position if and only if the current position is not a duplicate.
            //      i.e.: Check if the previous line had the same coordinates, and do NOT insert keys 
            //            of the previous line if it is the case.
            if self.insert_previous_keys()? {     
                continue 'line
            }
            // ---- Save the current coordinate for the next iteration if this position is not a duplicate.
            self.save_previous_line();
    

            // ---- Go to INFO field.
            self.skip_fields(5)?;                                      // 6 
            let info = self.reader.next_field()?.split(';').collect::<Vec<&str>>();

            // ---- Check if Bi-Allelic and skip line if not.
            if info.iter().any(|&field| field == "MULTI_ALLELIC") {
                self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
                self.coordinate_buffer.clear();
                continue 'line
            }

            // ---- Get the mutation type from INFO...
            let vtype = info.iter()
                .find(|&&field| field.starts_with("VT=")).unwrap()
                .split('=')
                .collect::<Vec<&str>>()[1];

            // ...and skip line if this is not an SNP.
            if vtype != "SNP" {
                self.clear_buffers()?;
                self.coordinate_buffer.clear();
                continue 'line
            }
            // ---- Extract population allele frequency.
            let mut pop_afs = info.iter()
                .filter(|&&field| field.contains("_AF"))
                .map(|field| field.split("_AF=").collect())
                .collect::<Vec<Vec<&str>>>();          
            pop_afs.sort();

            // ---- and insert the relevant entries within our frequency_set.
            for pop_af in pop_afs.iter() {
                let pop_af_key = format!("{}{} {}", std::str::from_utf8(&self.coordinate_buffer)?, pop_af[0], pop_af[1]);
                self.frequency_keys.push(pop_af_key);
            }


            self.skip_fields(1)?;                 // Skip INFO field
            self.fill_genotypes_buffer()?;        
            self.print_progress(50000)?;              
            self.parse_keys()?;
            self.genotypes_buffer.clear();
            self.coordinate_buffer.clear();
        }

        self.insert_previous_keys()?;
        Ok(())
    }

    /// Finish the construction of our genotype/frequency sets and flush their underlying writers.
    pub fn finish_build(self) -> Result<(), Box<dyn Error>> {
        self.builder.finish()?;
        self.frq_builder.finish()?;
        Ok(())
    }

    /// Parse the current VCF line's coordinate and fill `self.coordinate_buffer` with it.
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

    /// remove the last character of self.coordinate_buffer (i.e. the vcf field separator) and add our own Field separator.
    /// # Arguments:
    /// - `field`: byte character value of the desired field separator.
     fn add_field_separator(&mut self, field: u8) -> Result<(), Box<dyn Error>> {
        self.coordinate_buffer.pop();
        self.coordinate_buffer.push(field);
        Ok(())
    }
    
    /// Check if the current vcf coordinate is the same as the previous line.
    fn duplicate_line(&self) -> bool {
        self.coordinate_buffer == self.previous_coordinate
    }

    /// Checks if the current VCF line coordinate is a duplicate of the previous.
    /// - Skip and clear buffers if this is the case...
    /// - Flush the contents of our buffers within their respective SetBuilders if it is not the case 
    fn insert_previous_keys(&mut self) -> Result<bool, Box<dyn Error>> {
        if self.duplicate_line() {
            self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
            self.coordinate_buffer.clear();
            Ok(true)
        } else {
            self.insert_keys()?;       // If the line is different, we can then insert the keys of the previous line.
            Ok(false)
        }
    }

    /// Reset `self.reader` to the next line and clear our buffers.
    fn clear_buffers(&mut self) -> Result<(), Box<dyn Error>> {
        self.reader.source.read_until(b'\n', &mut Vec::new())?;
        self.genotype_keys.clear();
        self.frequency_keys.clear();
        Ok(())
    }

    /// Flush the contents of `self.genotypes_keys` and `self.frequency_keys` into their respective setbuilders.
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

    /// Clone and keep track of the current coordinate into `self.previous_coordinate`
    fn save_previous_line(&mut self) {
        self.previous_coordinate = self.coordinate_buffer.clone();
    }

    /// Skip a defined number of column fields within the VCF.
    /// # Arguments:
    /// - `n` number of fields to skip.
    fn skip_fields(&mut self, n: usize) -> Result<(), Box<dyn Error>> {
        for _ in 0..n {
            self.reader.source.read_until(b'\t', &mut Vec::new())?;
        }    
        Ok(())
    }

    /// read the contents of the whole vcf line and dump them into `self.genotypes_buffer`
    fn fill_genotypes_buffer(&mut self) -> Result<(), Box<dyn Error>> {   
        self.reader.source.read_until(b'\n', &mut self.genotypes_buffer)?;
        Ok(())
    }

    /// Log this Indexer's progress into the console, if `self.counter` is a multiple of `every_n`
    fn print_progress(&mut self, every_n: usize) -> Result<(), Box<dyn Error>> {
        if self.counter % every_n == 0 {
            info!("{: >9} {:?}", self.counter, std::str::from_utf8(&self.coordinate_buffer)?);
        }
        self.counter+= 1;
        Ok(())
    }

    /// Parse the genotypes of all the requested SampleTags (`self.samples.keys()`) and index them within `self.genotype_keys`
    unsafe fn parse_keys(&mut self) -> Result<(), Box<dyn Error>> {
        for sample_tag in self.samples.keys() {            

            // ---- Complete value to insert.
            let mut key = self.coordinate_buffer.clone();
            key.append(sample_tag.id().clone().as_mut_vec());  // UNSAFE: add sampleID
            key.push(b' ');

            // ---- Extract genotype info for thesample and add it to the key.
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

    // -------------------- Fetch the input vcf paths within `vcf_dir`
    info!("Fetching input VCF files in {}", &fst_cli.vcf_dir.to_str().unwrap());
    let mut input_vcf_paths = vcf::get_input_vcfs(&fst_cli.vcf_dir).unwrap();
    input_vcf_paths.sort();

    // --------------------- Fetch the input samples definition file
    let input_panel_path = match fst_cli.panel.clone() {
        Some(path) => path,
        None => vcf::fetch_input_panel(&fst_cli.vcf_dir)?,
    };

    // --------------------- Open the first vcf and assign column field indices for each sample
    //                       @ TODO: This could cause bugs if samples are not sorted in the same 
    //                               manner across multiple vcf files
    let mut panel = VCFPanelReader::new(&input_panel_path)?;
    panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;

    // --------------------- Subset the panel according to user-provided (super-)population id(s)
    match &fst_cli.pop_subset {
        Some(subset) => panel.subset_panel(subset),
        None         => (),
    }

    // --------------------- Setup a ThreadPool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(fst_cli.threads)
        .build()
        .unwrap();

    let output_panel_path = format!("{}/{}",
        fst_cli.output_dir.to_str().unwrap(),
        input_panel_path.file_name().unwrap().to_str().unwrap());
    info!("output_panel_path: {output_panel_path}");
    panel.copy_from_source(Path::new(&output_panel_path))?;


    pool.scope(|scope|{
        for vcf in input_vcf_paths.iter(){
            scope.spawn(|_| {
                // -------------------------- Format output file stem
                let file_stem = vcf.as_path()
                    .file_name()
                    .unwrap()
                    .to_str().unwrap()
                    .replace(".vcf.gz", "")
                    .replace(".vcf", "");

                // -------------------------- add population id to the output filename if the user provided some.
                let pop_tag = match &fst_cli.pop_subset {
                    Some(subset) => format!("-{}", subset.join("-")),
                    None => "".to_string(),
                };

                // -------------------------- Format the output file and destination directory.
                let output_path =format!("{}/{}{}", fst_cli.output_dir.to_str().unwrap(), file_stem, pop_tag);
                info!("Output path: {output_path:?}");

                // -------------------------- Build an FST set for this vcf file. 
                let mut setbuilder = VCFIndexer::new(vcf, &output_path, panel.into_transposed_btreemap(), fst_cli.decompression_threads).unwrap();
                unsafe {setbuilder.build_fst().unwrap();}

                // -------------------------- Finish the construction of our `.fst` and `.fst.frq` sets.
                setbuilder.finish_build().unwrap();
            });
        }
    });
    Ok(())
}
