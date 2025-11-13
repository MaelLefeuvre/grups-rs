use crate::snp::{Allele, ParseAlleleError};
use anyhow::Result;
use std::error::Error;

pub mod error;
pub use error::NucleotideError;

pub mod phred;
pub use phred::{Phred, PHRED_ASCII_BASE};


/// Simple struct representing a given nucleotides.
/// - base  : the nucleotide character -> generally preformatted by Pileup
/// - phred : Base-quality. Expressed in phred-33 scale.
/// 
/// # TODO: migrate nucleotide formating from `Pileup::new()` to `Nucleotide::new()`
#[derive(Debug, Clone, Copy)]
pub struct Nucleotide {
    pub base: Allele,
    pub phred: Phred,
}

impl Nucleotide {
    /// Instantiate a new `Nucleotide` struct, from a parsed pileup `base` and associated PHRED-33 `score` 
    #[must_use]
    pub fn new(base: Allele, score: impl Into<Phred>) -> Nucleotide {
        Nucleotide { base, phred: score.into() }
    }

    pub fn try_new<T>(base: T, score: char) -> Result<Nucleotide, NucleotideError> 
    where   T       : TryInto<Allele, Error = ParseAlleleError>,
            T::Error: Error + Send + Sync + 'static
    {
        let base = base.try_into().map_err(NucleotideError::ParseAllele)?;
        Ok(Nucleotide::new(base, score))
    }

    /// Convert the BQ score back to the ASCII format
    #[must_use]
    pub fn get_score_ascii(&self) -> char {
        self.phred.into()
    }

    /// Convert phred score to sequencing error probability
    #[must_use]
    pub fn error_prob(&self) -> f64 {
        // P = 10^(-Q/10)
        self.phred.as_prob()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn score_convert() -> Result<()> {
        for score in 0..PHRED_ASCII_BASE {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;
            let nucleotide = Nucleotide::try_new('A', ascii_score)?;
            assert_eq!(nucleotide.phred, Phred::new(ascii_score));
        }
        Ok(())
    }

    #[test]
    fn ascii_score_getter() -> Result<()> {
        for score in 0..PHRED_ASCII_BASE {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;
            let nucleotide = Nucleotide::try_new('A', ascii_score)?;
            assert_eq!(nucleotide.get_score_ascii(), ascii_score);
        }
        Ok(())
    }

    #[test]
    fn compute_error_prob() -> Result<()> {
        #![allow(clippy::float_cmp)]
        let scores         = [0,   10,  20,   30,    40    ];
        let expected_probs = [1.0, 0.1, 0.01, 0.001, 0.0001];

        for (i, score) in scores.iter().enumerate() {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;

            let nucleotide = Nucleotide::try_new('A', ascii_score)?;
            println!("{}", nucleotide.error_prob());
            assert_eq!(nucleotide.error_prob(), expected_probs[i]);
        }
        Ok(())
    }
}