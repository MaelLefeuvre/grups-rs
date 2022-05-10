use super::vcf::SampleTag;
use std::error::Error;

pub trait GenotypeReader {
    fn get_alleles(&self, sample_tag: &SampleTag )                  -> Option<[u8; 2]>;
    fn compute_local_cont_af(&self, contam_ind_ids: &[&SampleTag]) -> Result<f64, Box<dyn Error>>;
    fn get_pop_allele_frequency(&self, pop: &String)               -> Result<f64, Box<dyn Error>>;
}