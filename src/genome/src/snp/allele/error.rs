use thiserror::Error;

#[derive(Error, Debug)]
#[error("Failed to parse allele into a valid character")]
pub struct ParseAlleleError;