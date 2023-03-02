use thiserror::Error;

use super::{ChrIdxError, ParsePositionError};

#[derive(Error, Debug)]
pub enum CoordinateError {
    #[error("Failed to parse Coordinate because of an Invalid ChrIdx value")]
    ParseChrIdx(#[from] ChrIdxError),

    #[error("Failed to parse Coordinate because of an invalid Position value")]
    ParsePosition(#[from] ParsePositionError),

    #[error("Failed to parse Coordinate: missing delimiter '{0}' in string")]
    MissingDelimiter(char),

    #[error("Failed to parse Coordinate from byte slice.")]
    ParseFromSlice(#[from] std::array::TryFromSliceError)
}