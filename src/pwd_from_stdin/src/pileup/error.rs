use std::error::Error;
use std::fmt;

/// `PileupError` struct for error propagation and chaining.
///  - `RefSkip`       -> when char ['>', '<'] are found within a pileup string.
///                       used in `Pileup::new()`
///  - `UnequalLength` -> Nucleotide and Quality string, should have an equal length
///                       after filtration -> `LengthError` is raised if this is not the
///                       case.
///  - `ParseLine`     -> Not yet Implemented. General error which is raised if a 
///                       character failed to parse.
#[derive(Debug)]
pub enum PileupError {
    RefSkip,
    UnequalLength,
    //ParseLine,
}

impl Error for PileupError {}

impl fmt::Display for PileupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::RefSkip       => write!(f, "Cannot handle reference skips within pileup. ('>' '<')"),
            Self::UnequalLength => write!(f, "Length of nucleotides and scores differ vector differ."),
            //Self::ParseLine     => write!(f, "Failed to parse pileup line"),
        }
    }
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
        assert!( !error.is_empty())
    }

    //#[test]
    //fn pileup_err_parseline_display() {
    //    let error = format!("{}", PileupError::ParseLine);
    //    assert!( !error.is_empty())
    //}

}