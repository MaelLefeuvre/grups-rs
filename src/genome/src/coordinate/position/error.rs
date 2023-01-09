use thiserror::Error;

#[derive(Error, Debug)]
#[error("Failed to parse Physical Position into a valid u32: {0}")]
pub struct ParsePositionError(#[from] std::num::ParseIntError);