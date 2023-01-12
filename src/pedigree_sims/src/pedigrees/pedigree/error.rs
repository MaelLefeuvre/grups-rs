use thiserror::Error;

#[derive(Error, Debug)]
pub enum PedigreeError {

    #[error(transparent)]
    ParsePedigree(#[from] std::io::Error),

    #[error("Unexpected pedigree parameter")]
    EmptyParam,

    #[error("Failed to compare alleles")]
    FailedAlleleComparison,

    #[error("Missing or invalid contaminant list in pedigree")]
    MissingContaminant,

    #[error("Failed to retrieve a valid random sample from our panel")]
    MissingSampleTag,

    #[error("Failed to assign allele in individual {0}")]
    FailedAlleleAssignment(String),
}