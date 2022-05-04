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

pub struct FSTReader {
    genotypes_set: Set<Vec<u8>>,
    frequency_set: Set<Vec<u8>>,
    genotypes: HashMap<String, [u8; 2]>,
    frequencies: HashMap<String, f64>,
}

impl FSTReader {
    pub fn new(path: &str) -> Self {
        let genotypes_set = Self::get_set(&format!("{path}"));
        let frequency_set = Self::get_set(&format!("{path}.frq"));
        Self{genotypes_set, frequency_set, genotypes: HashMap::new(), frequencies: HashMap::new()}
    }

    fn get_set(path: &str) -> Set<Vec<u8>> {
        println!("{path}");
        let mut file_handle = File::open(path).unwrap();
        let mut bytes = vec![];
        std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
        Set::new(bytes).unwrap()
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

    pub fn get_pop_allele_frequency(&self, pop: &String) -> Option<&f64> {
        self.frequencies.get(pop)
    }

    pub fn get_alleles(&self, sample_id: &String ) -> Option<[u8; 2]> {
        Some(self.genotypes[sample_id].map(|all| all - 48))
    }

    pub fn has_genotypes(&self) -> bool {
       ! self.genotypes.is_empty()
    }

    pub fn compute_local_cont_af(&self, contam_ind_ids: &Vec<&SampleTag>) -> Result<f64, Box<dyn Error>>{
        let mut ref_allele_count = 0;
        let mut alt_allele_count = 0;
        for id in contam_ind_ids.iter().map(|tag| tag.id()) {
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

}