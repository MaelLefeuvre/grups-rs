
pub const PHRED_ASCII_BASE: u8 = 33;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Phred(u8);

impl From<char> for Phred {
    fn from(value: char) -> Self {
        Phred( (value as u8) - PHRED_ASCII_BASE )
    }
}

impl From<Phred> for char {
    fn from(val: Phred) -> Self {
        (val.0 + PHRED_ASCII_BASE) as char 
    }
}


impl From<u8> for Phred {
    fn from(value: u8) -> Self {
        Self(value)
    }
}

impl PartialEq<char> for Phred {
    fn eq(&self, other: &char) -> bool {
        self.0 == Phred::from(*other).0
    }
}

impl Phred {
    pub fn new(val: impl Into<Phred>) -> Self {
        val.into()
    }

    #[must_use]
    pub fn as_prob(&self) -> f64 {
        f64::powf(10.0, -f64::from(self.0) / 10.0) // P = 10^(-Q/10)
    }

    #[must_use]
    pub fn score(&self) -> u8 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_convert() {
        for raw_score in 0..PHRED_ASCII_BASE {
            let ascii_score = (raw_score + PHRED_ASCII_BASE) as char;
            let phred = Phred::new(ascii_score);
            assert_eq!(phred, Phred(raw_score));
        }
    }

    #[test]
    fn phred_from_char() {
        let raw = '@';
        let phred = Phred::from(raw);
        assert_eq!(phred.0, raw as u8 - PHRED_ASCII_BASE);
    }

    #[test]
    fn char_into_phred() {
        let raw = '@';
        let phred = Phred::new(raw);
        assert_eq!(raw, phred.into());
    }

    #[test]
    fn prob() {
        #![allow(clippy::float_cmp)]
        assert_eq!(Phred::from('+').as_prob(), 0.1);    // Phred 10 -> 10%
        assert_eq!(Phred::from('5').as_prob(), 0.01);   // PHRED 20 -> 1%
        assert_eq!(Phred::from('?').as_prob(), 0.001);  // PHRED 30 -> 0.1%
        assert_eq!(Phred::from('I').as_prob(), 0.0001); // PHRED 40 -> 0.01%
    }

    #[test] 
    fn char_equality() {
        for score in 0..PHRED_ASCII_BASE {
            assert_eq!(Phred::from(score), (score + PHRED_ASCII_BASE) as char);
        }
    }
}