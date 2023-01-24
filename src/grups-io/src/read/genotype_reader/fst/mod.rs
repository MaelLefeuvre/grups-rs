use std::{fs::File, path::{Path, PathBuf}, io::Read};

use crate::{
    parse,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader, GenotypeReaderError},
};


use genome::coordinate::Coordinate;
use located_error::prelude::*;

use memmap::Mmap;
use ahash::AHashMap;
use log::{info, debug};
use fst::{Set, IntoStreamer, automaton::{Automaton, Str}};

mod error;
use error::FSTReaderError;

pub const FST_EXT: &str = "fst";
pub const FRQ_EXT: &str = "fst.frq";

/// GenotypeReader from a pair of genotypes (`.fst`) and frequencies (`.fst.frq`) FST-Sets
/// #### Implemented Traits: GenotypeReader
pub struct FSTReader {
    genotypes_set: Set<Vec<u8>>,
    frequency_set: Set<Vec<u8>>,
    genotypes    : AHashMap<String, [u8; 2]>,
    frequencies  : AHashMap<String, f64>,
}


impl GenotypeReader for FSTReader {
    // Return the alleles for a given SampleTag. Search is performed using `sample_tag.id()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2]> {
        use GenotypeReaderError::MissingAlleles;
        let loc_msg = || format!("While retrieving alleles of {}", sample_tag.id());
        match self.genotypes.get(sample_tag.id()) {
            Some(alleles) => Ok(alleles.map(|all| all -48)),
            None          => Err(MissingAlleles).with_loc(loc_msg)
        }
    }

    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64> {
        use GenotypeReaderError::MissingFreq;
        match self.frequencies.get(pop) {
            Some(freq) => Ok(*freq),
            None       => Err(MissingFreq(pop.to_string())).loc("while parsing coordinate")
        }
    }

    fn fetch_input_files(input_dir: &Path) -> Result<Vec<PathBuf>> {
        let mut fsts = parse::fetch_input_files(input_dir, &[FST_EXT]).loc("While searching for candidate .fst files.")?;
        fsts.dedup();
        debug!("Found the following fst file candidates: {:#?}", fsts);
        // Ensure there's a matching 'fst.frq' file for each found '.fst' file
        for fst in fsts.iter() {
            Self::match_input_frq(fst).with_loc(|| format!("While fetching for .fst files in {}", input_dir.display()))?;
        }
        Ok(fsts)
    }

}

impl FSTReader {
    /// Instantiate a new reader from a `.fst` file.
    /// The companion `.fst.frq` file must be located within the same directory and display the same file-stem.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    pub fn new(path: &str) -> Result<Self> {
        let loc_msg = || format!("While attempting to create FSTReader from {path}");
        let genotypes_set = Self::get_set_memory(path).with_loc(loc_msg)?;
        let frequency_set = Self::get_set_memory(&format!("{path}.frq")).with_loc(loc_msg)?; //@TODO: const FRQ_EXT should be used in place.
        Ok(Self{genotypes_set, frequency_set, genotypes: AHashMap::new(), frequencies: AHashMap::new()})
    }

    /// Load a raw fst-set into memory and store it as a raw vector of bytes.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    #[allow(dead_code)]
    fn get_set_memory(path: &str) -> Result<Set<Vec<u8>>> {
        let loc_msg = || format!("While attempting to set the memory of {path}");
        use FSTReaderError::{OpenFST, LoadFST, BuildFST};
        info!("Loading in memory : {path}");
        // ---- Get file size
        let size = std::fs::metadata(path).map_err(OpenFST).with_loc(loc_msg)?.len();
        // ---- Open FST file in ro mode.
        let mut file_handle = File::open(path).map_err(OpenFST).with_loc(loc_msg)?;

        // ---- Load FST set into memory.
        let mut bytes = Vec::with_capacity(size as usize);
        Read::read_to_end(&mut file_handle, &mut bytes).map_err(LoadFST).with_loc(loc_msg)?;

        // ---- Generate FST set
        Set::new(bytes).map_err(BuildFST).with_loc(loc_msg)
    }

    /// Load a raw fst-set as a memory-mapped file.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    #[allow(dead_code)]
    fn get_set_memmap(path: &str) -> Set<Mmap> {
        info!("Loading mem-map: {path}");
        let mmap = unsafe { Mmap::map(&File::open(path).unwrap()).unwrap() };
        Set::new(mmap).unwrap()
    }

    /// Public wrapper for `find_chromosome()`. returns a sorted, unduplicated list of chromosomes contained within the set.
    pub fn find_chromosomes(&self) -> Result<Vec<u8>> {
        let root = self.genotypes_set.as_fst().root();
        let mut string = Vec::with_capacity(4);
        let mut chromosomes = Vec::new();
        self.find_chromosome(root, &mut string, &mut chromosomes)
            .loc("Failed to determine which chromosome(s) were contained within the set. File is possibly corrupt.")?;
        chromosomes.sort_unstable();
        chromosomes.dedup();
        Ok(chromosomes)
    }


    /// Recursively search through the raw genotype set and populate `chromosomes` each time a chromosome has been found.
    /// # Arguments:
    /// - `node`       : a raw fst-index node.
    /// - `string`     : raw byte-string. This vector is passed and constructed upon each recursion step.
    /// - `chromosomes`: vector of chromosomes names (in u8) form. This is the final desired output.
    ///                  Each entry of `chromosomes` is constructed from `string`.
    fn find_chromosome(&self, node: fst::raw::Node, string: &mut [u8], chromosomes: &mut Vec<u8>) -> Result<()> {
        use FSTReaderError::{InvalidUTF8, InvalidChrId};
        let loc_msg = "While recursing through the available chromosomes within the FST set";
        let sep = b' '; // FST-set fields are space-separated. 
        for transition in node.transitions() {
            if transition.inp != sep { // If the next character is not a space, add the current character to the string being constructed.
                let mut local_string = string.to_owned();
                local_string.push(transition.inp); 
                self.find_chromosome(self.genotypes_set.as_fst().node(transition.addr), &mut local_string, chromosomes)?;
            }
            else { // If we've found a separator character, the current string is fully defined. -> add this string to our chromosome list.
                
                // ---- Parse the string into a valid chromosome u8 value.
                let chromosome = std::str::from_utf8(string)
                    .map_err(InvalidUTF8)
                    .and_then(|c| c.parse::<u8>().map_err(InvalidChrId))
                    .loc(loc_msg)?;

                chromosomes.push(chromosome);
            }
        }
        Ok(())

    }

    /// Recursively search through the index for a given chromosome and return true if it has been found.
    /// # Arguments:
    /// - `chr`: chromosome name to look for.
    pub fn contains_chr(&self, chr: u8) -> bool {
        let regex= format!("{chr} ");
        let mut node = self.genotypes_set.as_fst().root();
        for b in regex.as_bytes() {
            match node.find_input(*b) {
                None => return false,
                Some(i) => {
                    node = self.genotypes_set.as_fst().node(node.transition_addr(i));
                }
            }
        }
        true
    }

    fn format_coordinate_pattern(coordinate: &Coordinate) -> String {
        format!("{} {:0>9}", coordinate.chromosome, coordinate.position)
    }
    /// Populate the `genotypes` field using genome coordinates.
    /// `self.genotypes` is thus a HashMap, with keys==<sample-id>, and values==<alleles>
    /// # Arguments:
    /// - `chr`: chromosome name
    /// - `pos`: 0-based chromosome position
    pub fn search_coordinate_genotypes(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        // ---- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let regex = Self::format_coordinate_pattern(coordinate);
        let matcher = Str::new(&regex).starts_with();

        // ---- Search through the set and format each match within `self.genotypes` (key=<sample-id>, val=<alleles>)
        self.genotypes_set.search(&matcher)
            .into_stream()
            .into_bytes()
            .iter()
            .map(|string| {
                string.split(|char| *char == b' ')
                .skip(2)
            })
            .for_each(|mut v| {
                let tag = unsafe { std::str::from_utf8_unchecked(v.next().unwrap())};
                let alleles: [u8; 2] = v.next()
                    .map(|slice| {
                        slice.try_into().unwrap()
                    }).unwrap();
                self.genotypes.insert(tag.to_string(), alleles);
            });
    }

    /// Populate the `frequencies` field using genome coordinates.
    /// `self.frequencies` is thus a HashMap, with keys==<sample-id>, and values==<allele-frequency>
    /// # Arguments:
    /// - `chr`: chromosome name
    /// - `pos`: 0-based chromosome position
    pub fn search_coordinate_frequencies(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        //----- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let regex = Self::format_coordinate_pattern(coordinate);
        let matcher = Str::new(&regex).starts_with();

        // ---- Search through the set and format each match within `self.frequencies` (key=<sample-id>, val=<allele-frequency>)
        self.frequency_set.search(&matcher)
            .into_stream()
            .into_bytes()
            .iter()
            .map(|string| {
                string.split(|char| *char == b' ')
                .skip(2)
            })
            .for_each(|mut v| {
                let pop = unsafe { std::str::from_utf8_unchecked(v.next().unwrap())};
                let freq: f64 = unsafe { 
                    std::str::from_utf8_unchecked(v.next().unwrap())
                    .parse::<f64>()
                    .unwrap()
                };
                self.frequencies.insert(pop.to_string(), freq);
            });
    }


    /// Clear `self.genotypes` and `self.frequencies` buffer.
    /// Generally called before skipping to the next coordinate.
    pub fn clear_buffers(&mut self) {
        self.genotypes.clear();
        self.frequencies.clear();
    }

    /// Check if the `self.genotypes` field has been filled.
    pub fn has_genotypes(&self) -> bool {
       ! self.genotypes.is_empty()
    }

    /// Search for a companion `.fst.frq` file next to the provided `.fst`. Returns an error if there are none.
    /// # Arguments: 
    /// - `fst_file`: path leading to a `.fst` file.
    /// # Errors:
    /// - returns an error if there are no `.fst.frq` companion file for the provided `fst_file`
    fn match_input_frq(fst: &Path) -> Result<()> {
        let frq = fst.with_extension(FRQ_EXT);
        if ! frq.exists() {
            return Err(FSTReaderError::MatchFrqFile{path: fst.to_path_buf()})
                .loc("Note that .fst files must have a matching .fst.frq file within the same location.")
        }
        Ok(())
    }
    
}