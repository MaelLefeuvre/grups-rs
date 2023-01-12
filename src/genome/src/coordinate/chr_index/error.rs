use thiserror::Error;

#[derive(Error, Debug)]
#[error("Failed to parse Chromosome name into a valid u8")]
pub struct ChrIdxError(#[from] std::num::ParseIntError);