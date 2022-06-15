use super::vcf::SampleTag;
use std::error::Error;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait GenotypeReader {
    /// Extract the alleles from the GenotypeReader, given a SampleTag
    fn get_alleles(&self, sample_tag: &SampleTag )                     -> Option<[u8; 2]>;
    /// Extract the population allele frequency of the current line from the GenotypeReader, given a population id.
    fn get_pop_allele_frequency(&self, pop: &str)                      -> Result<f64, Box<dyn Error>>;
}