use thiserror::Error;

#[derive(Error, Debug)]
pub enum NucleotideError {
    #[error("Failed to parse nucleotide base into a valid Allele")]
    ParseAllele(#[from] crate::snp::ParseAlleleError)
}