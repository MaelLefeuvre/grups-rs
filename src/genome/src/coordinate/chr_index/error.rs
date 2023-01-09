use thiserror::Error;

#[derive(Error, Debug)]
#[error("Failed to parse Chromosome name into a valid u8: {0}")]
pub struct ChrIdxError(pub String);