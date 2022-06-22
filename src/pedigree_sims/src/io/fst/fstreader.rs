use std::{
    fs::File,
    collections::HashMap,
    error::Error
};

use fst::{
    Set,
    IntoStreamer,
    automaton::{Automaton, Str},
};

use crate::{io::vcf::SampleTag};
use crate::io::genotype_reader::GenotypeReader;

use log::info;
use memmap::Mmap;

/// GenotypeReader from a pair of genotypes (`.fst`) and frequencies (`.fst.frq`) FST-Sets
/// #### Implemented Traits: GenotypeReader
pub struct FSTReader {
    genotypes_set: Set<Vec<u8>>,
    frequency_set: Set<Vec<u8>>,
    genotypes: HashMap<String, [u8; 2]>,
    frequencies: HashMap<String, f64>,
}

impl FSTReader {

    /// Instantiate a new reader from a `.fst` file.
    /// The companion `.fst.frq` file must be located within the same directory and display the same file-stem.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    pub fn new(path: &str) -> Self {
        let genotypes_set = Self::get_set_memory(path);
        let frequency_set = Self::get_set_memory(&format!("{path}.frq"));
        Self{genotypes_set, frequency_set, genotypes: HashMap::new(), frequencies: HashMap::new()}
    }

    /// Load a raw fst-set into memory and store it as a raw vector of bytes.
    /// # Arguments:
    /// - `path`: path leading to the targeted `.fst` file.
    #[allow(dead_code)]
    fn get_set_memory(path: &str) -> Set<Vec<u8>> {
        info!("Loading in memory : {path}");
        let mut file_handle = File::open(path).unwrap();
        let mut bytes = vec![];
        std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
        Set::new(bytes).unwrap()
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

    /// public wrapper for `find_chromosome()`. returns the list of chromosomes contained within the set.
    pub fn find_chromosomes(&self) -> Vec<u8> {
        let root = self.genotypes_set.as_fst().root();
        let mut string = Vec::new();
        let mut chromosomes = Vec::new();
        self.find_chromosome(root, &mut string, &mut chromosomes);
        chromosomes
    }


    /// Recursively search through the raw genotype set and populate `chromosomes` each time a chromosome has been found.
    /// # Arguments:
    /// - `node`       : a raw fst-index node.
    /// - `string`     : raw byte-string. This vector is passed and constructed upon each recursion step.
    /// - `chromosomes`: vector of chromosomes names (in u8) form. This is the final desired output.
    ///                  Each entry of `chromosomes` is constructed from `string`.
    fn find_chromosome(&self, node: fst::raw::Node, string: &mut [u8], chromosomes: &mut Vec<u8>) {
        let sep = b' '; // FST-set fields are space-separated. 
        for transition in node.transitions() {
            if transition.inp != sep { // If the next character is not a space, add the current character to the string being constructed.
                let mut local_string = string.to_owned();
                local_string.push(transition.inp); 
                self.find_chromosome(self.genotypes_set.as_fst().node(transition.addr), &mut local_string, chromosomes);
            }
            else { // If we've found a separator character, the current string is fully defined. -> add this string to our chromosome list.
                let chromosome = std::str::from_utf8(string).unwrap().parse::<u8>().unwrap();
                chromosomes.push(chromosome);
            }
        }

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

    /// Populate the `genotypes` field using genome coordinates.
    /// `self.genotypes` is thus a HashMap, with keys==<sample-id>, and values==<alleles>
    /// # Arguments:
    /// - `chr`: chromosome name
    /// - `pos`: 0-based chromosome position
    pub fn search_coordinate_genotypes(&mut self, chr: u8, pos: u32) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        //----- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let regex = format!("{chr} {pos:0>9}");
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
    pub fn search_coordinate_frequencies(&mut self, chr: u8, pos: u32) {
        // ---- Create a new fst-matcher, matching any entry starting with our chromosome coordinates.
        //----- genotype-fst index fields are '{chr} {pos} {id} {alleles}'
        let regex = format!("{chr} {pos:0>9}");
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
}

impl GenotypeReader for FSTReader {
    // Return the alleles for a given SampleTag. Search is performed using `sample_tag.id()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Option<[u8; 2]> {
        Some(self.genotypes[sample_tag.id()].map(|all| all - 48))
    }

    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, Box<dyn Error>> {
        match self.frequencies.get(pop) {
            Some(freq) => Ok(*freq),
            None             => Err(format!("Missing allele frequency in .frq file for pop {} at this coordinate.", pop).into())
        }
    }
}