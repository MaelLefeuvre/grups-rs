use std::{
    ops::Deref
};

#[derive(Debug, Clone)]
pub struct Allele {
    pos    : u32,
    allele : u8,
    af     : f64,
}

impl Allele {
    pub fn new(pos: u32, allele: u8, af: f64) -> Self {
        Allele{pos, allele, af}
    }
    pub fn get_pos(&self) -> u32 {self.pos}

    pub fn get_allele(&self) -> u8 {self.allele}

    pub fn get_af(&self) -> f64 {self.af}
}

#[derive(Debug, Clone)]
pub struct Alleles (Vec<Allele>);

impl Alleles{
    pub fn add_allele(&mut self, allele: Allele){
        self.0.push(allele);
    }
}

impl Default for Alleles {
    fn default() -> Self {
        Alleles(Vec::new())
    }
}

impl Extend<Allele> for Alleles {
    fn extend<T: IntoIterator<Item=Allele>>(&mut self, iter: T) {
        for elem in iter {
            self.add_allele(elem);
        }
    }
}

impl Deref for Alleles {
    type Target = Vec<Allele>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
