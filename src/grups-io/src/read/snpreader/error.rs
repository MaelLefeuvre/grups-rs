use thiserror::Error;
use genome::coordinate::Coordinate;

use super::SNPREADER_VALID_FILE_FORMATS;
fn common_help_msg() -> String {
    format!("Please provide '--targets' with either one of the accepted file formats: {SNPREADER_VALID_FILE_FORMATS:?}")
}
#[derive(Error, Debug)]
pub enum SNPReaderError {
    #[error("Cannot handle target file format: {0}, {}", common_help_msg())]
    InvalidFileFormat(String),

    #[error("The provided targets file is missing a file extension. {}", common_help_msg())]
    MissingExtension,

    #[error("Failed to open file {0}")]
    OpenFile(String, #[source] std::io::Error),

    #[error("Target SNP coordinate {0} is missing its REF and/or ALT allele, when it is required.")]
    MissingAltRef(Coordinate),

    #[error(
        "Cannot filter-out transitions if reference and/or alternate alleles is unknown. \
        Please provide '--targets' with a file containing known reference and alternate alleles."
    )]
    UnknownAlleles,
}