use thiserror::Error;

/// `PileupError` struct for error propagation and chaining.
///  - `RefSkip`:  when char `['>', '<']` are found within a pileup string. Used in `Pileup::new()`
/// 
///  - `UnequalLength`:  Nucleotide and Quality string, should have an equal length after filtration.
///    `LengthError` is raised if this is not the case.
/// 
///  - `ParseLine`: Not yet Implemented. General error which is raised if a character failed to parse.
/// 
#[derive(Error, Debug)]
pub enum PileupError {
    #[error("Cannot handle reference skips within pileup. ('>' '<')")]
    RefSkip,
    #[error("Length of nucleotides and scores differ vector differ.")]
    UnequalLength,

    #[error("Failed to parse allele within pileup")]
    ParseRef(#[from] genome::snp::ParseAlleleError),

    #[error("Failed to parse chromosome into a valid u8")]
    ParseChr(#[from] genome::coordinate::ChrIdxError),
    
    #[error("Failed to parse genomic position into a valid u32")]
    ParsePos(#[from] genome::coordinate::ParsePositionError),

    #[error("Failed to parse pileup depth int a valid u8")]
    ParseDepth(#[from] std::num::ParseIntError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pileup_err_refskip_display() {
        let error = format!("{}", PileupError::RefSkip);
        assert!( !error.is_empty());
    }

    #[test]
    fn pileup_err_unequal_display() {
        let error = format!("{}", PileupError::UnequalLength);
        assert!( !error.is_empty());
    }
}