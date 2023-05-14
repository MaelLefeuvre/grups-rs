use thiserror::Error;
use genome::coordinate::Coordinate;

fn locate_genomic_position(coord: &[u8]) -> String {
    match Coordinate::try_from(&coord[0..5]).ok() {
        Some(coord) => format!("{coord}"),
        None        => "[None]".to_string()
    }
}


#[derive(Error, Debug)]
pub enum GenomeFstError {
    #[error("Failed to open target output file")]
    CreateFile(#[source] std::io::Error),

    #[error("Failed to instantiate FST Set Builder")]
    CreateSetBuilder(#[source] fst::Error),

    #[error("The provided file or directory cannot be parsed into a string")]
    DisplayPath,

    #[error("Failed to instantiate threadpool")]
    BuildThreadPool(#[source] rayon::ThreadPoolBuildError),

    #[error("{} Missing VT Info Tag", locate_genomic_position(c))]
    MissingVTTag{c: Vec<u8>},

    #[error("{} Missing EOL in entry", locate_genomic_position(c))]
    MissingEOL{c: Vec<u8>},

    #[error("{} Failed to parse population allele frequency to f32", locate_genomic_position(c))]
    ParseAlleleFrequency{c: Vec<u8>},

    #[error("{} INFO field appears to contain duplicate population allele frequency tags", locate_genomic_position(c))]
    DuplicatePopFreqTag{c: Vec<u8>},

    #[error("{} Failed to parse or encode chromosome", locate_genomic_position(c))]
    EncodeChr{c: Vec<u8>},

    #[error("{} Failed to parse or encode position", locate_genomic_position(c))]
    EncodePos{c: Vec<u8>},

    #[error("{} Failed to fill buffer with the next field.", locate_genomic_position(c))]
    ReadField{c: Vec<u8>},

    #[error("{} Failed to fill genotypes buffer.", locate_genomic_position(c))]
    FillGenotypes{c: Vec<u8>},

    #[error("Failed to insert previous entry in the FST set.")]
    InsertPreviousEntry,

    #[error("{} Failed to insert key within the Fst", locate_genomic_position(c))]
    InsertFstKey{c: Vec<u8>},

    #[error("Failed to finish FST set build")]
    CompleteBuild(#[source] fst::Error),
}