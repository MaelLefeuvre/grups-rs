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

pub struct FSTReader {
    genotypes_set: Set<Vec<u8>>,
    frequency_set: Set<Vec<u8>>,
    genotypes: HashMap<String, [u8; 2]>,
    frequencies: HashMap<String, f64>,
}

impl FSTReader {
    pub fn new(path: &str) -> Self {
        let genotypes_set = Self::get_set_memory(path);
        let frequency_set = Self::get_set_memory(&format!("{path}.frq"));
        Self{genotypes_set, frequency_set, genotypes: HashMap::new(), frequencies: HashMap::new()}
    }

    #[allow(dead_code)]
    fn get_set_memory(path: &str) -> Set<Vec<u8>> {
        info!("Loading in memory : {path}");
        let mut file_handle = File::open(path).unwrap();
        let mut bytes = vec![];
        std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
        Set::new(bytes).unwrap()
    }

    #[allow(dead_code)]
    fn get_set_memmap(path: &str) -> Set<Mmap> {
        info!("Loading mem-map: {path}");
        let mmap = unsafe { Mmap::map(&File::open(path).unwrap()).unwrap() };
        Set::new(mmap).unwrap()
    }

    pub fn find_chromosomes(&self) -> Vec<u8> {
        let node = self.genotypes_set.as_fst().root();
        let mut string = Vec::new();
        let mut chromosomes = Vec::new();
        self.find_chromosome(node, &mut string, &mut chromosomes);
        chromosomes
    }


    pub fn find_chromosome(&self, node: fst::raw::Node, string: &mut Vec<u8>, chromosomes: &mut Vec<u8>) {
        let sep = b' ';
        match node.find_input(sep) {
            None => {
                for transition in node.transitions() {
                    string.push(transition.inp);
                    let mut local_string = string.clone();
                    self.find_chromosome(self.genotypes_set.as_fst().node(transition.addr), &mut local_string, chromosomes);
                }

            },
            Some(_) => {
                let chromosome = std::str::from_utf8(&string).unwrap().parse::<u8>().unwrap();
                chromosomes.push(chromosome);

            },
        }
    }

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

    pub fn search_coordinate_genotypes(&mut self, chr: u8, pos: u32) {
        let regex = format!("{chr} {pos:0>9}");
        let matcher = Str::new(&regex).starts_with();
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

    pub fn search_coordinate_frequencies(&mut self, chr: u8, pos: u32) {
        let regex = format!("{chr} {pos:0>9}");
        let matcher = Str::new(&regex).starts_with();
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


    pub fn clear_buffers(&mut self) {
        self.genotypes.clear();
        self.frequencies.clear();
    }

    pub fn has_genotypes(&self) -> bool {
       ! self.genotypes.is_empty()
    }
}

impl GenotypeReader for FSTReader {
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Option<[u8; 2]> {
        Some(self.genotypes[sample_tag.id()].map(|all| all - 48))
    }

    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, Box<dyn Error>> {
        match self.frequencies.get(pop) {
            Some(freq) => Ok(*freq),
            None             => Err(format!("Missing allele frequency in .frq file for pop {} at this coordinate.", pop).into())
        }
    }
}