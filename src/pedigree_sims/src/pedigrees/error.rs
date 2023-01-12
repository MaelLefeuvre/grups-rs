use thiserror::Error;

#[derive(Debug, Error)]
pub enum PedigreeError {
    #[error("Failed to access and update previous position tracker.")]
    InvalidCoordinate,   

    #[error("Attempting to mutably access a missing pedigree vector using the comparison label '{0}' as a key.")]
    MissingPedVec(String),

    #[error("Pedigree Vector does not contain any contaminant")]
    MissingContaminant
}