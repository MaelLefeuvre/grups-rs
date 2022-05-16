const PHRED_ASCII_BASE: u8 = 33;

/// Simple struct representing a given nucleotides.
/// - base  : the nucleotide character -> generally preformatted by Pileup
/// - phred : Base-quality. Expressed in phred-33 scale.
/// 
/// # TODO: migrate nucleotide formating from `Pileup::new()` to `Nucleotide::new()`
#[derive(Debug)]
pub struct Nucleotide {
    pub base: char,
    pub phred: u8,
}

impl Nucleotide {
    pub fn new(base: char, score: char) -> Nucleotide {
        let phred = Nucleotide::to_phred(score);       
        Nucleotide {base: base, phred}
    }

    /// Convert the BQ score back to the ASCII format
    pub fn get_score_ascii(&self) -> char {
        (self.phred + PHRED_ASCII_BASE) as char 
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

    #[test]
    fn score_convert(){
        for score in 0..PHRED_ASCII_BASE {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;
            let nucleotide = Nucleotide::new('A', ascii_score);
            assert_eq!(nucleotide.phred, score);
        }
    }

    #[test]
    fn ascii_score_getter(){
        for score in 0..PHRED_ASCII_BASE {
            let ascii_score = (score + PHRED_ASCII_BASE) as char;
            let nucleotide = Nucleotide::new('A', ascii_score);
            assert_eq!(nucleotide.get_score_ascii(), ascii_score);
        }
    }
}