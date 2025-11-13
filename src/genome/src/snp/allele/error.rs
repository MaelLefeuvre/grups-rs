use thiserror::Error;

#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseAlleleError {

    #[error("Failed to parse character into a valid Allele")]
    InvalidCharacter,
    
    #[error("Failed to parse IUPAC Ambiguity code into a valid Allele.")]
    IUPACAmbiguityCodeCharacterFound
}