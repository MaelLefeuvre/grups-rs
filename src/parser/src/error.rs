use thiserror::Error;
use crate::FileEntity;

#[derive(Error, Debug)]
pub enum ParserError{
    #[error("Invalid slice or value format for --{arg}. [{err}]")]
    ParseArg{arg: String, err: String},

    #[error("{0}")]
    ParseRange(&'static str),

    #[error("--min-depth must be greater than 0")]
    InsufficientDepthError,

    #[error("Neither --pileup, nor the stdin buffer are being sollicited.")]
    MissingPileupInput,

    #[error("{0} {1} does not exist")]
    MissingFileEntity(FileEntity, String),

    #[error("{1} is not a {0}")]
    InvalidFileEntity(FileEntity, String),

    #[error("The provided value must lie between {0} and {1}")]
    ParseRatio(f64, f64),

    #[error("Failed to generate an output file prefix. Note that file prefixes are generated from the input pileup filestem")]
    ParseOutputPrefix,

    #[error("{0} already exists. Use --overwrite to force.")]
    CannotOverwrite(String)

}