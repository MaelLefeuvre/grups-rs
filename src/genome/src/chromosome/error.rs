use thiserror::Error;
#[derive(Error, Debug)]
pub enum ChromosomeError {
    #[error("Failed to parse chromosome length")]
    ParseLength,

    #[error("Failed to parse chromosome name")]
    ParseChrIdx
}


