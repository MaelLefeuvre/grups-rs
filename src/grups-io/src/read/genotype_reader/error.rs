use thiserror::Error;

#[derive(Error, Debug)]
pub enum GenotypeReaderError {
    #[error("Failed to retrieve the alleles of the individual")]
    MissingAlleles,

    #[error("Failed to retrieve the allele frequency of population {0}")]
    MissingFreq(String),

    #[error("Invali")]
    InvalidSampleIndex
}