use thiserror::Error;

#[derive(Debug, Error)]
pub enum IndividualError {
    #[error("Strand index is unexpectedly empty.")]
    MissingStrands,

    #[error("Attempting to access empty alleles.")]
    MissingAlleles,

    #[error("Failed to assign alleles of individual")]
    InvalidAlleleAssignment,

    #[error("Failed to assign sex of individual")]
    InvalidSexAssignment,

    #[error("Invalid allele assignment: Assigning strands is bound to fail, as this individual has no parents.")]
    MissingParents,

    #[error("Unknown or Missing Sex in individual")]
    UnknownOrMissingSex,

    #[error("Invalid or Spurious recombination event")]
    InvalidOrSpuriousRecombinationEvent,

    #[error("Spurious allele assignment: {alleles:?}")]
    SpuriousAlleleAssignment{alleles: [u8; 2]}
}