use std::path::PathBuf;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ParseError {
    #[error("Failed to create parent directory")]
    CreateParentDirectory(#[source] std::io::Error),

    #[error("File or directory returned an empty string, and may contain invalid UTF-8 characters")]
    InvalidFilename,

    #[error("'{}' already exists within the output directory. Use '--overwrite' to force, or specify a different output directory with '--output-dir'", path.display())]
    OverwriteDisallowed{path: PathBuf},

    #[error("Could not find any valid input file within the requested directory: {}", dir.display())]
    MissingInput{dir: PathBuf}
}