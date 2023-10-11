use genome::coordinate::Coordinate;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ComparisonError {
    #[error("Cannot find corresponding Jackknife Block for the provided coordinate {0}")]
    MissingBlock(Coordinate),

    #[error("Invalid Pileup index '{0}' for individual '{1}'. The number of individuals within your pileup file may be lower than expected...
    Ensure that each specified sample does not exceed the number of individuals within the pileup file.
    Also note that indices specified with --sample are zero-based indices."
    )]
    InvalidPileupIndex(usize, String),
}