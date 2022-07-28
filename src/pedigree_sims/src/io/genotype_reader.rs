use super::vcf::SampleTag;

#[cfg(test)]
use mockall::{automock, predicate::*};

#[cfg_attr(test, automock)]
pub trait GenotypeReader {
    /// Extract the alleles from the GenotypeReader, given a SampleTag
    fn get_alleles(&self, sample_tag: &SampleTag )                     -> Result<[u8; 2], String>;
    /// Extract the population allele frequency of the current line from the GenotypeReader, given a population id.
    fn get_pop_allele_frequency(&self, pop: &str)                      -> Result<f64, String>;
}