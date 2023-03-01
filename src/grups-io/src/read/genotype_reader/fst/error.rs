use std::path::PathBuf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum FSTReaderError {
    #[error("Failed to open the source FST file.")]
    OpenFST(#[source] std::io::Error),

    #[error("Failed to read the contents of the source FST file.")]
    LoadFST(#[source] std::io::Error),

    #[error("Failed to construct a valid set from the source FST file (data might be corrupted)")]
    BuildFST(#[source] fst::Error),

    #[error("Could not find any valid paired '.fst.frq' file for {}", path.display())]
    MatchFrqFile{path: PathBuf}

}