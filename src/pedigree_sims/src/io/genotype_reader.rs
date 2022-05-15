use crate::contaminant::Contaminant;

use super::vcf::SampleTag;
use std::error::Error;

pub trait GenotypeReader {
    fn get_alleles(&self, sample_tag: &SampleTag )                     -> Option<[u8; 2]>;
    fn compute_local_cont_af(&self, contaminant: Option<&Contaminant>) -> Result<f64, Box<dyn Error>>;
    fn get_pop_allele_frequency(&self, pop: &str)                      -> Result<f64, Box<dyn Error>>;
}