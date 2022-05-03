use std::{
    fs::File,
    collections::HashMap,
};

use fst::{
    Set,
    IntoStreamer,
    automaton::{Automaton, Str},
};

pub struct FSTReader {
    pub genotypes_set: Set<Vec<u8>>,
    pub frequency_set: Set<Vec<u8>>
}

impl FSTReader {
    pub fn new(path: &str) -> Self {
        let genotypes_set = Self::get_set(&format!("{path}"));
        let frequency_set = Self::get_set(&format!("{path}.frq"));
        Self{genotypes_set, frequency_set}
    }

    fn get_set(path: &str) -> Set<Vec<u8>> {
        println!("{path}");
        let mut file_handle = File::open(path).unwrap();
        let mut bytes = vec![];
        std::io::Read::read_to_end(&mut file_handle, &mut bytes).unwrap();
        Set::new(bytes).unwrap()
    }

    pub fn search_coordinate_genotypes<'a>(&'a self, chr: u8, pos: u32) -> HashMap<String, [u8; 2]> {
        let mut output = HashMap::new();
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
                    
                output.insert(tag, alleles);
            });
        output
    }

    pub fn search_coordinate_frequencies(&self, chr: u8, pos: u32)  -> HashMap<String, f64> {
        let mut output = HashMap::new();
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
                output.insert(pop, alleles);
            });
        output
    }
}