use genome::coordinate::Coordinate;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ComparisonError {
    #[error("Cannot find corresponding Jackknife Block for the provided coordinate {0}")]
    MissingBlock(Coordinate),
}