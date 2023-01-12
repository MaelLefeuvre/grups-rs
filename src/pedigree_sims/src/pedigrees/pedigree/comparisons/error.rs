use thiserror::Error;

#[derive(Error, Debug)]
pub enum ComparisonError {
    #[error("Failed to compare alleles of comparison.")]
    CompareAllele,

    #[error("Failed to sample a random allele during simulations")]
    SampleAllele,

    #[error("Failed to select a random erroneous base when simulating sequencing error.")]
    SimSeqError,
}