use std::{fs::File, path::{Path, PathBuf}, io::Read};

use crate::{
    parse,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader, GenotypeReaderError},
};


use genome::coordinate::Coordinate;
use located_error::prelude::*;

use memmap2::Mmap;
use ahash::AHashMap;
use log::{info, debug};
use fst::{Set, IntoStreamer, automaton::{Automaton, Str, StartsWith}, Streamer};

mod error;
use error::FSTReaderError;

pub const FST_EXT: &str = "fst";
pub const FRQ_EXT: &str = "fst.frq";

/// GenotypeReader from a pair of genotypes (`.fst`) and frequencies (`.fst.frq`) FST-Sets
/// #### Implemented Traits: GenotypeReader
pub struct FSTReader {
    genotypes_set: Set<Vec<u8>>,
    frequency_set: Set<Vec<u8>>,
    genotypes    : AHashMap<u128, [u8; 2]>,
    frequencies  : AHashMap<String, f32>,
}


impl GenotypeReader for FSTReader {
    // Return the alleles for a given SampleTag. Search is performed using `sample_tag.id()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2]> {
        use GenotypeReaderError::MissingAlleles;
        let loc_msg = || format!("While retrieving alleles of {}", sample_tag.id());

        match self.genotypes.get(&sample_tag.hashed_id()) {
            Some(alleles) => Ok(alleles.map(|all| all -48)),
            None          => Err(MissingAlleles).with_loc(loc_msg)
        }
    }

    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f32> {
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
        debug!("File size (bytes): {size}");
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
        let mut chromosomes = Vec::from_iter( root.transitions().map(|transition| transition.inp) );
        chromosomes.sort_unstable();
        chromosomes.dedup();
        Ok(chromosomes)
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

    /// @TODO: This should be a macro ?
    #[inline]
    fn format_coordinate_pattern(coord_bytes: &[u8; 5]) -> StartsWith<Str> {
        let regex = unsafe { std::str::from_utf8_unchecked(coord_bytes) };
        Str::new(regex).starts_with()
    }
    
    /// Populate the `genotypes` field using genome coordinates.
    /// `self.genotypes` is thus a HashMap, with keys==<sample-id>, and values==<alleles>
    #[inline]
    pub fn search_coordinate_genotypes(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        // ---- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let coord_bytes: [u8; 5] = (*coordinate).into();
        let matcher = Self::format_coordinate_pattern(&coord_bytes);

        // ---- Search through the set and format each match within `self.genotypes` (key=<sample-id>, val=<alleles>)
        let mut stream = self.genotypes_set.search(&matcher).into_stream();
        while let Some(key) = stream.next() {
            // ---- Retrieve alleles and sample id.
            let key = &mut key[5..].to_vec();              // Skip coordinate bytes.
            let (id, alleles) = key.split_at(key.len()-2); // alleles are the last two u8.

            // ---- Hash sample id. to u128
            let hashed_tag = SampleTag::hash_id_u128(id);
            
            let alleles: [u8; 2] = alleles.try_into().unwrap();
            self.genotypes.insert(hashed_tag, alleles);
        }
    }

    /// Populate the `frequencies` field using genome coordinates.
    /// `self.frequencies` is thus a HashMap, with keys==<sample-id>, and values==<allele-frequency>
    /// # Arguments:
    /// - `chr`: chromosome name
    /// - `pos`: 0-based chromosome position
    #[inline]
    pub fn search_coordinate_frequencies(&mut self, coordinate: &Coordinate) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        //----- genotype-fst index fields are '{chr(u8)}{pos(u32_be)}{pop(chars)}{freq(f32_be)}'
        let coord_bytes: [u8; 5] = (*coordinate).into();
        let matcher = Self::format_coordinate_pattern(&coord_bytes);

        // ---- Search through the set and format each match within `self.frequencies` (key=<sample-id>, val=<allele-frequency>)
        let mut stream = self.frequency_set.search(&matcher).into_stream();
        while let Some(key) = stream.next() {
            // ---- Retrieve population tag and frequency.
            let key              = &mut key[5..].to_vec(); // Go from the end of the coordinate to the end.
            let (tag, freq_be)   = key.split_at(key.len()-4); 

            // ---- Parse population tag
            let pop              = unsafe {String::from_utf8_unchecked(tag.to_vec())};

            // ---- Convert frequency to f32
            let freq_be: [u8; 4] = freq_be.try_into().unwrap();
            let freq: f32        = f32::from_be_bytes(freq_be);

            self.frequencies.insert(pop, freq);
        }
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