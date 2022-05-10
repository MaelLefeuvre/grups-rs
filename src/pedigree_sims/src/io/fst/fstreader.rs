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

use crate::io::vcf::SampleTag;
use crate::io::genotype_reader::GenotypeReader;

use log::info;
use memmap::Mmap;

pub struct FSTReader {
    genotypes_set: Set<Mmap>,
    frequency_set: Set<Mmap>,
    genotypes: HashMap<String, [u8; 2]>,
    frequencies: HashMap<String, f64>,
}

impl FSTReader {
    pub fn new(path: &str) -> Self {
        let genotypes_set = Self::get_set_memmap(path);
        let frequency_set = Self::get_set_memmap(&format!("{path}.frq"));
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

    fn get_set_memmap(path: &str) -> Set<Mmap> {
        info!("Loading mem-map: {path}");
        let mmap = unsafe { Mmap::map(&File::open(path).unwrap()).unwrap() };
        Set::new(mmap).unwrap()
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
            .into_strs().unwrap()
            .iter()
            .map(|string| string.split(' ').skip(2).collect::<Vec<&str>>())
            .for_each(|v| {
                let tag: String = v[0].to_string();
                let alleles: [u8; 2] = v[1].chars().map(|c| c as u8)
                    .collect::<Vec<u8>>()
                    .try_into()
                    .unwrap_or_else(|v: Vec<u8>| panic!("Expected a Vec of length {} but it was {}", 2, v.len()));
                    
                self.genotypes.insert(tag, alleles);
            });
    }

    pub fn search_coordinate_frequencies(&mut self, chr: u8, pos: u32) {
        let regex = format!("{chr} {pos:0>9}");
        let matcher = Str::new(&regex).starts_with();
        self.frequency_set.search(&matcher)
            .into_stream()
            .into_strs().unwrap()
            .iter()
            .map(|string| string.split(' ').skip(2).collect::<Vec<&str>>())
            .for_each(|v| {
                let pop: String = v[0].to_string();
                let alleles: f64 = v[1].parse().unwrap();
                self.frequencies.insert(pop, alleles);
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

    fn compute_local_cont_af(&self, contam_ind_ids: &[&SampleTag]) -> Result<f64, Box<dyn Error>> {
        let mut ref_allele_count = 0;
        let mut alt_allele_count = 0;
        for id in contam_ind_ids.iter() {
            let cont_alleles = self.get_alleles(id).unwrap();
            for allele in cont_alleles.into_iter() {
                match allele {
                    0 => ref_allele_count +=1,
                    1 => alt_allele_count +=1,
                    _ => panic!("Contaminating individual is multiallelic: {allele}")
                }
            }
        }
        Ok((alt_allele_count as f64) /(alt_allele_count as f64 + ref_allele_count as f64))
    }

    fn get_pop_allele_frequency(&self, pop: &String) -> Result<f64, Box<dyn Error>> {
        match self.frequencies.get(pop) {
            Some(freq) => Ok(*freq),
            None             => Err(format!("Missing allele frequency in .frq file for pop {} at this coordinate.", pop).into())
        }
    }
}