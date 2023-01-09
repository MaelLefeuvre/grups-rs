use genome::snp::Allele;
use anyhow::Result;
use std::error::Error;

const PHRED_ASCII_BASE: u8 = 33;

/// Simple struct representing a given nucleotides.
/// - base  : the nucleotide character -> generally preformatted by Pileup
/// - phred : Base-quality. Expressed in phred-33 scale.
/// 
/// # TODO: migrate nucleotide formating from `Pileup::new()` to `Nucleotide::new()`
#[derive(Debug, Clone, Copy)]
pub struct Nucleotide {
    pub base: Allele,
    pub phred: u8,
}

impl Nucleotide {
    /// Instantiate a new `Nucleotide` struct, from a parsed pileup `base` and associated PHRED-33 `score` 
    #[must_use]
    pub fn new(base: Allele, score: char) -> Nucleotide {
        let phred = Nucleotide::to_phred(score);       
        Nucleotide {base, phred}
    }

    pub fn try_new<T>(base: T, score: char) -> Result<Nucleotide> 
    where   T: TryInto<Allele>,
            T::Error: Error + Send + Sync + 'static
    {
        let phred = Nucleotide::to_phred(score);
        let base = base.try_into()?;
        Ok(Nucleotide{base, phred})
    }

    /// Convert the BQ score back to the ASCII format
    #[must_use]
    pub fn get_score_ascii(&self) -> char {
        (self.phred + PHRED_ASCII_BASE) as char 
    }

    /// Convert phred score to sequencing error probability
    #[must_use]
    pub fn error_prob(&self) -> f64 {
        // P = 10^(-Q/10)
        f64::powf(10.0, -1.0 * f64::from(self.phred) / 10.0)
    }

    /// Convert an ASCII BQ score to phred-33.
    /// Used during construction.
    fn to_phred(score: char) -> u8 {
        (score as u8) - PHRED_ASCII_BASE
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
            assert_eq!(nucleotide.phred, score);
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
        let scores          = [0,   10,  20,   30,    40    ];
        let expected_probs = [1.0, 0.1, 0.01, 0.001, 0.0001];

        for (i, score) in scores.iter().enumerate() {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;

            let nucleotide = Nucleotide::try_new('A', ascii_score)?;
            println!("{}", nucleotide.error_prob());
            assert_eq!(nucleotide.error_prob(), expected_probs[i])
        }
        Ok(())
    }
}