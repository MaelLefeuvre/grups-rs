use thiserror::Error;

#[derive(Error, Debug)]
pub enum InfoFieldError {
    #[error("INFO field does not contain any 'VT=' tag.")]
    MissingVT,

    #[error("INFO field does not contain any '{0}' tag.")]
    MissingAF(String),

    #[error("Failed to parse population allele frequency to a float")]
    ParseAlleleFrequency(#[source] std::num::ParseFloatError),
}