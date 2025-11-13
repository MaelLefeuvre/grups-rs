use std::{collections::BTreeMap, fs::File, io::{BufRead, BufWriter}, path::Path, process};

use genome::coordinate::Coordinate;

use located_error::prelude::*;
use parser::VCFFst;
use grups_io::read::{
    SampleTag, PanelReader,
    genotype_reader::{GenotypeReader, VCFReader, FST_EXT, FRQ_EXT},
};

use fst::SetBuilder;
use log::{info, error, trace, debug};
use rayon::{self}; 

mod error;
use error::GenomeFstError;


pub enum AFStrategy {
    ComputeFromAlleles,
    ParseFromVcf
}


///Create a new `SetBuilder` from a file
fn create_set_builder(file: &str) -> Result<SetBuilder<BufWriter<File>>> {
    use GenomeFstError::{CreateFile, CreateSetBuilder};
    let loc_msg = "While attempting to generate a new Set Builder";
    let file    = File::create(file).map_err(CreateFile).loc(loc_msg)?;
    SetBuilder::new(BufWriter::new(file)).map_err(CreateSetBuilder).loc(loc_msg)
}

/// Convert a `.vcf(.gz)` file into a Finite-State-Transducer Set.
/// See the following resources for a recap on the theory and applications of FST Sets :
/// - `<https://cs.nyu.edu/~mohri/pub/fla.pdf` (DOI: 10.1007/978-3-540-39886-8_29) 
/// - `<https://blog.burntsushi.net/transducers/>`
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
/// - `builder`            : FST `SetBuilder` for genotypes (this is what constructs our `.fst` file)
/// - `frq_builder`        : FST `SetBuilder` for allele frequencies (this is what constructs ou `.fst.frq` file)
/// - `samples`            : Transposed input panel definition map, sorted across `SampleTag.id()`
/// - `counter`            : Keeps track of the number of processed VCF lines.
/// - `coordinate_buffer`  : Coordinate of the current VCF line
/// - `previous_coordinate`: Coordinate of the previous VCF line
/// - `genotypes_buffer`   : Genotypes of the current VCF line
/// - `genotypes_keys`     : Temporary buffer of the values that should be inserted into the genotypes set builder
/// - `frequency_keys`     : Temporary buffer of the values that should be inserted into the frequency set builder
pub struct VCFIndexer<'a, 'panel>{
    reader              : VCFReader<'a>,
    gen_builder         : SetBuilder<BufWriter<File>>,
    frq_builder         : SetBuilder<BufWriter<File>>,
    samples             : BTreeMap<&'panel SampleTag, Vec<&'panel str>>,
    counter             : usize,
    coordinate_buffer   : Vec<u8>,
    previous_coordinate : Vec<u8>,
    genotypes_buffer    : Vec<u8>,
    genotype_keys       : Vec<Vec<u8>>,
    frequency_keys      : Vec<Vec<u8>>
}      

impl<'a, 'panel> VCFIndexer<'a, 'panel> {
    /// Initialize a new `VCFIndexer`.
    /// # Arguments:
    /// - `vcf`    : path leading to the target `.vcf`|`.vcf.gz` file
    /// - `samples`: transposed input samples definition file (sorted across `SampleTag.id()``)
    /// - `threads`: number of requested decompression threads (only relevant in the case of BGZF compressed `.vcf.gz` files)
    pub  fn new(vcf: &'a Path, output_file: &str, samples: BTreeMap<&'panel SampleTag, Vec<&'panel str>>, threads: usize) -> Result<Self> {
        let loc_msg = |ext| format!("While attempting to create a new VCFIndexer for {output_file}.{ext}");
        // ----------------------- Initialize Genotypes and Frequency Writers
        let gen_builder = create_set_builder(&format!("{output_file}.{FST_EXT}")).with_loc(|| loc_msg(FST_EXT))?;
        let frq_builder = create_set_builder(&format!("{output_file}.{FRQ_EXT}")).with_loc(|| loc_msg(FRQ_EXT))?;

        // ----------------------- Initialize VCFReader.
        let reader  = VCFReader::new(vcf, threads)?;
        Ok(Self {
            reader,
            gen_builder,
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
    pub fn build_fst(&mut self, allele_frequency_strategy: &AFStrategy) -> Result<()> {
        use GenomeFstError::MissingVTTag;
        let loc_msg = "While constructing FST set.";
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
            let info = self.reader.next_field()?.split(';').map(ToString::to_string).collect::<Vec<String>>();

            // ---- Check if Bi-Allelic and skip line if not.
            if info.iter().any(|field| field == "MULTI_ALLELIC") {
                self.clear_buffers()?;     // If so, skip this lines, and don't insert the keys of the previous line.
                self.coordinate_buffer.clear();
                continue 'line
            }

            // ---- Get the mutation type from INFO...
            let vtype = info.iter()
                .find(|&field| field.starts_with("VT="))
                .with_loc(||MissingVTTag {c:self.coordinate_buffer.clone()})?
                .split('=')
                .collect::<Vec<&str>>()[1];

            // ...and skip line if this is not an SNP.
            if vtype != "SNP" {
                self.clear_buffers().with_loc(||"While attempting to skip non-SNP position.")?;
                self.coordinate_buffer.clear();
                continue 'line
            }


            // ---- fill our genotypes buffer with the sample's alleles.
            self.skip_fields(1).with_loc(|| loc_msg)?;                 // Skip INFO field
            self.fill_genotypes_buffer().with_loc(|| loc_msg)?;        
            
            // ---- Extract or compute population allele frequency from INFO.
            let pop_afs = match allele_frequency_strategy {
                AFStrategy::ComputeFromAlleles => self.compute_allele_frequencies()?,
                AFStrategy::ParseFromVcf       => self.parse_allele_frequencies(&info)?,
            };
            
            debug!("{} | Pop_afs: {pop_afs:?}", Coordinate::try_from(&self.coordinate_buffer[0..5])?);

           // ---- and insert the relevant entries within our frequency_set.
            self.insert_frequency_keys(&pop_afs);
            
            self.print_progress(50000).with_loc(|| loc_msg)?;              
            self.parse_keys().with_loc(|| loc_msg)?;
            
            self.genotypes_buffer.clear();
            self.coordinate_buffer.clear();
        }

        self.insert_previous_keys()?;
        Ok(())
    }

    /// Finish the construction of our genotype/frequency sets and flush their underlying writers.
    /// # Errors
    /// - if either the `.fst` or `.fst.frq` file fails to finish building and writing itself.
    pub fn finish_build(self) -> Result<()> {
        use GenomeFstError::CompleteBuild;
        self.gen_builder.finish().map_err(CompleteBuild)?;
        self.frq_builder.finish().map_err(CompleteBuild)?;
        Ok(())
    }

    /// Parse the current VCF line's coordinate and fill `self.coordinate_buffer` with it.
    #[inline]
    fn fill_current_position_buffer(&mut self) -> Result<()> {
        use GenomeFstError::{ReadField, EncodeChr, EncodePos};
        // ---------------------------- Get current chromosome and position.
        let mut chr_buffer = Vec::with_capacity(3);
        self.reader.source.read_until(b'\t', &mut chr_buffer).with_loc(||ReadField { c: self.coordinate_buffer.clone() })?;
        chr_buffer.pop();

        let mut pos_buffer = Vec::with_capacity(10);
        self.reader.source.read_until(b'\t', &mut pos_buffer).with_loc(||ReadField { c: self.coordinate_buffer.clone() })?;
        pos_buffer.pop();

        let chr: u8 = str::from_utf8(&chr_buffer).ok()
            .and_then(|c| {
                match c.parse::<u8>() {
                    Ok(b) => Some(b),
                    Err(_) => {
                        match c == "X" || c == "chrX" {
                            true =>   Some(b'X'),
                            false =>  None
                        }
                    }
                }
            }).with_loc(||EncodeChr {c: self.coordinate_buffer.clone() })?;

        let pos: [u8; 4] = str::from_utf8(&pos_buffer).ok().and_then(|c| c.parse::<u32>().ok()).map(u32::to_be_bytes)
            .with_loc(||EncodePos { c: self.coordinate_buffer.clone() })?;

        self.coordinate_buffer.push(chr);
        self.coordinate_buffer.extend(pos);
        Ok(())
    }

    
    /// Check if the current vcf coordinate is the same as the previous line.
    #[inline]
    fn duplicate_line(&self) -> bool {
        self.coordinate_buffer == self.previous_coordinate
    }

    /// Checks if the current VCF line coordinate is a duplicate of the previous.
    /// - Skip and clear buffers if this is the case...
    /// - Flush the contents of our buffers within their respective `SetBuilders` if it is not the case 
    #[inline]
    fn insert_previous_keys(&mut self) -> Result<bool> {
        use GenomeFstError::InsertPreviousEntry;
        if self.duplicate_line() {
            // If so, skip this lines, and don't insert the keys of the previous line.
            self.clear_buffers().with_loc(||InsertPreviousEntry)?;     
            self.coordinate_buffer.clear();
            Ok(true)
        } else {
             // If the line is different, we can then insert the keys of the previous line.
            self.insert_keys().with_loc(||InsertPreviousEntry)?;      
            Ok(false)
        }
    }

    #[inline]
    fn insert_frequency_keys(&mut self, pop_afs: &BTreeMap<&str, f32>) {
        // ---- and insert the relevant entries within our frequency_set.
        for (pop_tag, pop_af) in pop_afs {
            // ---- Parse to: "{chromosome(u7)}{position(u32_be)}{pop(chars)}{freq(f32_be)}"
            let pop_tag    = pop_tag.as_bytes();
            let pop_af_be  = pop_af.to_be_bytes();
    
            let frequency_key = self.coordinate_buffer.iter()
                .chain(pop_tag.iter())
                .chain(pop_af_be.iter())
                .copied();
            self.frequency_keys.push(frequency_key.collect::<Vec<_>>());
        }
    }

    /// Reset `self.reader` to the next line and clear our buffers.
    #[inline]
    fn clear_buffers(&mut self) -> Result<()> {
        use GenomeFstError::MissingEOL;
        self.reader.source.read_until(b'\n', &mut Vec::new())
            .with_loc(||MissingEOL{c: self.coordinate_buffer.clone()})?;
        self.genotype_keys.clear();
        self.frequency_keys.clear();
        Ok(())
    }

    /// Flush the contents of `self.genotypes_keys` and `self.frequency_keys` into their respective setbuilders.
    #[inline]
    fn insert_keys(&mut self) -> Result<()> {
        use GenomeFstError::InsertFstKey;
        for gen in &self.genotype_keys {
            self.gen_builder.insert(gen).with_loc(|| InsertFstKey{c:self.coordinate_buffer.clone()})?;
        }
        self.genotype_keys.clear();

        for frq in &self.frequency_keys {
            self.frq_builder.insert(frq).with_loc(|| InsertFstKey{c:self.coordinate_buffer.clone()})?;
        } 
        self.frequency_keys.clear();
        Ok(())
    }

    /// Clone and keep track of the current coordinate into `self.previous_coordinate`
    #[inline]
    fn save_previous_line(&mut self) {
        self.previous_coordinate.clone_from(&self.coordinate_buffer);
    }

    /// Skip a defined number of column fields within the VCF.
    /// # Arguments:
    /// - `n` number of fields to skip.
    #[inline]
    fn skip_fields(&mut self, n: usize) -> Result<()> {
        use GenomeFstError::ReadField;
        for _ in 0..n {
            self.reader.source
                .read_until(b'\t', &mut Vec::new())
                .with_loc(||ReadField { c: self.coordinate_buffer.clone() })?;
        }    
        Ok(())
    }

    /// read the contents of the whole vcf line and dump them into `self.genotypes_buffer`
    #[inline]
    fn fill_genotypes_buffer(&mut self) -> Result<usize> {
        use GenomeFstError::FillGenotypes;
        self.reader.source.read_until(b'\n', &mut self.genotypes_buffer)
            .with_loc(||FillGenotypes{ c: self.coordinate_buffer.clone() })
    }

    /// Attempt to parse the current coordinate buffer as a valid coordinate
    #[inline]
    fn current_coordinate(&self) -> Result<Coordinate> {
        Ok(Coordinate::try_from(&self.coordinate_buffer[0..5]).map_err(GenomeFstError::InvalidCoordinate)?)
    }

    /// Log this Indexer's progress into the console, if `self.counter` is a multiple of `every_n`
    #[inline]
    fn print_progress(&mut self, every_n: usize) -> Result<()> {
        if self.counter % every_n == 0 {
            let coord = self.current_coordinate().loc("While attempting to print progress")?;
            info!("{: >9} {coord}", self.counter);
        }
        self.counter+= 1;
        Ok(())
    }

    /// Parse the genotypes of all the requested `SampleTags` (`self.samples.keys()`) and index them within `self.genotype_keys`
    #[inline]
    fn parse_keys(&mut self) -> Result<()>{

        let sample_genotypes: Vec<&[u8]> = self.genotypes_buffer.split(|c| *c == b'\t' || *c == b'\n').collect();
        for sample_tag in self.samples.keys() {            

            // ---- Complete value to insert.
            let mut key = self.coordinate_buffer.clone(); // @TODO: That clone can be removed?
            unsafe { key.append(sample_tag.id().clone().as_mut_vec()) };  // UNSAFE: add sampleID

            // ---- Extract genotype info for the sample and add it to the key.
            let geno_idx = sample_tag.idx().expect("Missing sample tag index"); 
            let sample_genotype = sample_genotypes.get(geno_idx)
                .with_loc(|| GenomeFstError::InvalidGenotypeIndex { 
                    c: self.coordinate_buffer[0..5].to_vec(), s: (sample_tag.id()).clone(), i: geno_idx
                })?;

            trace!("  - {} | {}:{}",
                self.current_coordinate().loc("While parsing keys")?,
                sample_tag.id(),
                str::from_utf8(sample_genotype).loc("While parsing keys")?
            );

            let (left, right) = match sample_genotype.len() {
                3 => Ok((sample_genotype[0], sample_genotype[2])), // Autosomal or Pseudo-autosomal region
                1 => Ok((sample_genotype[0], sample_genotype[0])), // X-chromosome (for male samples)
                o => Err(GenomeFstError::InvalidGenotypeKey{c: self.coordinate_buffer[0..5].to_vec(), s: (sample_tag.id()).clone(), l: o})
            }?;

            key.extend_from_slice(&[left, right]); // haplo1

            self.genotype_keys.push(key);
        }
        Ok(())
    }

    #[inline]
    fn parse_allele_frequencies<'info >(&self, info: &'info [String]) -> Result<BTreeMap<&'info str, f32>> {
        use GenomeFstError::{ParseAlleleFrequency, DuplicatePopFreqTag};
        let mut frequencies : BTreeMap<&str, f32> = BTreeMap::new();
        
        let pop_afs: Option<Vec<(&str, &str)>> = info.iter()
            .filter(|&field| field.contains("_AF"))
            .map(|field| field.split_once("_AF="))
            .collect();

        let pop_afs = pop_afs.with_loc(|| ParseAlleleFrequency{c:self.coordinate_buffer.clone()})?;
        for (pop, allele_frequency) in &pop_afs {
            if frequencies.contains_key(pop) {
                return Err(DuplicatePopFreqTag{c: self.coordinate_buffer.clone()})
                    .loc("While attempting to parse allele frequencies from INFO field")

            }
            frequencies.insert(pop, allele_frequency.parse::<f32>().map_err(|_| GenomeFstError::DisplayPath)?);
        } 

        Ok(frequencies)
    }

    #[inline]
    fn compute_allele_frequencies(&self) -> Result<BTreeMap<&'panel str, f32>> {
        let mut frequencies : BTreeMap<&str, f32> = BTreeMap::new();
        let mut allele_count: BTreeMap<&str, f32> = BTreeMap::new();

        let sample_genotypes: Vec<&[u8]> = self.genotypes_buffer.split(|c| *c == b'\t' || *c == b'\n').collect();

        for (sample_tag, pop_tags) in &self.samples {
            let genotype_index = sample_tag.idx().expect("Missing sample tag index");
            let sample_genotype = sample_genotypes.get(genotype_index)
                .with_loc(|| GenomeFstError::InvalidGenotypeIndex { 
                    c: self.coordinate_buffer[0..5].to_vec(), s: (sample_tag.id()).clone(), i: genotype_index
                })?;
            
            let (left, right) = match sample_genotype.len() {
                3  => Ok((sample_genotype[0], sample_genotype[2])), // Autosomal or Pseudo-autosomal region
                1  => Ok((sample_genotype[0], sample_genotype[0])), // X-chromosome (for male samples)
                o  => Err(GenomeFstError::InvalidGenotypeKey{c: self.coordinate_buffer[0..5].to_vec(), s: (sample_tag.id()).clone(), l: o})
            }?;

            if sample_genotypes.len() > 1 {
                for pop_tag in pop_tags { // Autosome, Pseudo-autosomal region, or female X-chromosome
                    *frequencies.entry(pop_tag).or_default()  += if left == 49 {1.0} else {0.0};
                    *frequencies.entry(pop_tag).or_default()  += if right == 49 {1.0} else {0.0};
                    *allele_count.entry(pop_tag).or_default() += 2.0;
                }
            } else {
                for pop_tag in pop_tags { // Male X-chromosome
                    *frequencies.entry(pop_tag).or_default()  += if left == 49 {1.0} else {0.0};
                    *allele_count.entry(pop_tag).or_default() += 1.0;
                }
            }

        }

        for (pop_tag, alternative_allele_count) in &mut frequencies {
            const FREQ_ROUND_DECIMAL: u32 = 4;
            let observed_alleles = allele_count.get(pop_tag).expect("Missing sample tag index");
            trace!("{pop_tag}: {alternative_allele_count} / {observed_alleles} = {:.4}", *alternative_allele_count / *observed_alleles);

            *alternative_allele_count /= observed_alleles;

            let round_factor = 10u32.pow(FREQ_ROUND_DECIMAL) as f32;
            *alternative_allele_count = (round_factor * *alternative_allele_count).round() / round_factor;
        }

        Ok(frequencies)
    }
}


pub fn run(fst_cli: &VCFFst) -> Result<()> {
    let loc_msg = "While attempting to initialize FST Set Builder";

    // -------------------- Fetch the input vcf paths within `vcf_dir`
    info!("Fetching input VCF files in {}", &fst_cli.vcf_dir.to_string_lossy());
    let mut input_vcf_paths = VCFReader::fetch_input_files(&fst_cli.vcf_dir).loc(loc_msg)?;
    input_vcf_paths.sort();

    // --------------------- Fetch the input samples definition file
    let mut panel = match fst_cli.panel.as_ref() {
        Some(path) => PanelReader::new(path),
        None       => PanelReader::from_dir(&fst_cli.vcf_dir),
    }.loc(loc_msg)?;
    // --------------------- Open the first vcf and assign column field indices for each sample
    //                       @ TODO: This could cause bugs if samples are not sorted in the same 
    //                               manner across multiple vcf files
    panel.assign_vcf_indexes(input_vcf_paths[0].as_path()).loc(loc_msg)?;

    // --------------------- Subset the panel according to user-provided (super-)population id(s)
    if let Some(subset) = &fst_cli.pop_subset { panel.subset_panel(subset) }

    // --------------------- Parse the output panel file.
    let Some(output_dir) = fst_cli.output_dir.to_str() else {
        return Err(GenomeFstError::DisplayPath).loc(loc_msg)
    };

    let Some(Some(source_file)) = panel.source_file.file_name().map(|path| path.to_str()) else {
        return Err(GenomeFstError::DisplayPath).loc(loc_msg)
    };
    let output_panel_path = format!("{output_dir}/{source_file}");
    
    info!("Output_panel_path: {output_panel_path}");
    panel.copy_from_source(Path::new(&output_panel_path))?;

    // --------------------- Check whether the user requested allele frequency recalculation or not.
    let frequency_strategy = match fst_cli.compute_pop_afs {
        true  => AFStrategy::ComputeFromAlleles,
        false => AFStrategy::ParseFromVcf,
    };

    // --------------------- Setup a ThreadPool
    info!("Setting up thread pool with {} threads", fst_cli.threads);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(fst_cli.threads)
        .build()
        .map_err(GenomeFstError::BuildThreadPool)
        .loc(loc_msg)?;

    pool.scope(|scope|{
        for vcf in &input_vcf_paths{
            scope.spawn(|_| {
                // -------------------------- Format output file stem
                let file_stem = vcf.as_path()
                    .file_name()
                    .and_then(|f| f.to_str())
                    .expect("Invalid vcf file name (non-UTF8 characters present?)")
                    .replace(".vcf.gz", "")
                    .replace(".vcf", "");

                // -------------------------- add population id to the output filename if the user provided some.
                let pop_tag = match &fst_cli.pop_subset {
                    Some(subset) => format!("-{}", subset.join("-")),
                    None => String::new(),
                };

                // -------------------------- Format the output file and destination directory.
                let output_path =format!("{output_dir}/{file_stem}{pop_tag}");
                info!("Output path: {output_path:?}");

                // -------------------------- Build an FST set for this vcf file. 
                let mut setbuilder = match VCFIndexer::new(
                        vcf,
                        &output_path,
                        panel.into_transposed_btreemap(),
                        fst_cli.decompression_threads
                    ).with_loc(|| "While constructing VCFIndexer") {
                        Ok(builder) => builder,
                        Err(e) => {
                            error!("{e:?}");
                            process::exit(1);
                        }
                    };

                if let Err(e) = setbuilder.build_fst(&frequency_strategy).with_loc(|| "While building FST from VCF") {
                    error!("{e:?}");
                    process::exit(1);
                }

                // -------------------------- Finish the construction of our `.fst` and `.fst.frq` sets.
                if let Err(e) = setbuilder.finish_build().with_loc(|| "While Finalizing FST build") {
                    error!("{e:?}");
                    process::exit(1);
                }
            });
        }
    });

    Ok(())
}
