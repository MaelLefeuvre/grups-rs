use crate::alleles::Alleles;

pub struct Chromatid {
    pub index   : usize, 
    pub name    : u8,
    pub length  : u32,
    pub alleles : Alleles
}

impl Chromatid {
    pub fn new(index: usize, name: u8, length: u32, alleles: Alleles) -> Self {
        Chromatid{index, name, length, alleles}
    }
}
