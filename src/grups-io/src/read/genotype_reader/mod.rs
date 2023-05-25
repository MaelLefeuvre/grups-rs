use std::path::{PathBuf, Path};

use crate::read::SampleTag;

mod vcf;
pub use vcf::VCFReader;

pub mod fst;
pub use self::fst::FSTReader;
pub use self::fst::{FST_EXT, FRQ_EXT};


mod error;
pub use error::GenotypeReaderError;
use anyhow::Result;

//#[cfg(test)]
use mockall::{automock, predicate::*};



// @ TODO: this mocking should really be a cfg(test), attribute, but since I'm 
// using it in pedigree_sims::pedigrees::contaminants, for testing purposes, I'm having trouble making MockGenotypeReader visible.
//#[cfg_attr(test, automock)]
#[automock]       
pub trait GenotypeReader {
    /// Extract the alleles from the GenotypeReader, given a SampleTag
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2]>;

    /// Extract the population allele frequency of the current line from the GenotypeReader, given a population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f32>;

    /// Iterate over the contetns of an OS-Directory and search for all the required input files.
    fn fetch_input_files(input_dir: &Path) -> Result<Vec<PathBuf>> where Self: Sized;
}
