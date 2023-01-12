use thiserror::Error;

#[derive(Debug, Error)]
pub enum IndividualError {
    #[error("Strand index is unexpectedly empty.")]
    MissingStrands,

    #[error("Attempting to access empty alleles.")]
    MissingAlleles,

    #[error("Failed to assign alleles of individual")]
    InvalidAlleleAssignment,

    #[error("Invalid allele assignment: Assigning strands is bound to fail, as this individual has no parents.")]
    MissingParents,
}