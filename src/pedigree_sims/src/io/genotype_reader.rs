use super::vcf::SampleTag;
use std::error::Error;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait GenotypeReader {
    fn get_alleles(&self, sample_tag: &SampleTag )                     -> Option<[u8; 2]>;
    fn get_pop_allele_frequency(&self, pop: &str)                      -> Result<f64, Box<dyn Error>>;
}