use std::path::PathBuf;
use thiserror::{Error, self};


#[derive(Error, Debug)]
pub enum GeneticMapError {
    #[error("Failed to parse '{}' into a valid genetic map", map.display())]
    ParseMap{map: PathBuf},

    #[error("Failed to read contents of the provided directory")]
    ReadDir,

    #[error("Failed to find or parse any genetic-map in the provided directory.")]
    EmptyDir,

    #[error("Line {0} appears to be invalid")]
    InvalidLine(usize),

    #[error("Failed to parse chromosome field @ line {0}")]
    ParseChr(usize),

    #[error("Failed to parse position field @ line {0}")]
    ParsePos(usize),

    #[error("Failed to parse genetic rate field @ line {0}")]
    ParseRate(usize),

    #[error("File appears to be missing a field @ line {0}")]
    InvalidFields(usize)
}
