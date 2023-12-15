use thiserror::Error;

use super::PedFormatField;

#[derive(Debug, Error)]
pub enum PedigreeBuilderError{
    #[error("Failed to open template pedigree definition file: {0}")]
    OpenFile(#[source] std::io::Error),

    #[error("Encountered IO error when reading line n°{lineno} of the pedigree definition file: {source}")]
    IoError{source: std::io::Error, lineno: usize},
    
    #[error("Invalid number of fields within the provided template pedigree file. Must be 3, 4 or 6. Got {0}")]
    InvalidFieldNumber(usize),

    #[error("Invalid header within the provided template pedigree definition file. Check the first line of this file.")]
    InvalidHeader,

    #[error("Failed to retrieve field n°{0}")]
    MissingField(PedFormatField),

    #[error("Failed to retrieve {0} field at index {1}")]
    RetrieveField(PedFormatField, usize),

    #[error("Failed to add individual {0} while parsing line n°{1} in the pedigree definition file")]
    AddIndividual(String, usize),

}
