mod error;
use std::{fmt::{Display, Debug}, str::FromStr, borrow::Borrow};

use error::ParseAlleleError;

/// Struct representing a pileup Allele
/// Usual cases:     A: Adenine    C: Cytosine    G: Guanine    T: Thymine
/// Specific cases:  N: Unknown    D: Deletion (CIGAR D)
#[derive(Debug, Hash, PartialEq, Eq, Clone, Copy)]
pub enum Allele { A, C, G, T, N, D}

impl From<&Allele> for char {
    fn from(value: &Allele) -> Self {
        match value {
            Allele::A   => 'A',
            Allele::C   => 'C',
            Allele::G   => 'G',
            Allele::T   => 'T',
            Allele::N   => 'N',
            Allele::D   => '*',
        }
    }
}

impl TryFrom<char> for Allele {
    type Error = ParseAlleleError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        use self::Allele::*;
        match value {
            'A'       => Ok(A),
            'C'       => Ok(C),
            'G'       => Ok(G),
            'T'       => Ok(T),
            'N'       => Ok(N),
            '*'       => Ok(D),
             _  => Err(ParseAlleleError)
        }
    }
}

impl FromStr for Allele {
    type Err = ParseAlleleError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let char = s.parse::<char>().map_err(|_| ParseAlleleError)?;
        Ok(Self::try_from(char)?)
    }
}

impl Display for Allele {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&char::from(self) , f)
    }
}

impl Borrow<char> for Allele {
    fn borrow(&self) -> &char {
        match self {
            Self::A => &'A',
            Self::C => &'C',
            Self::G => &'G',
            Self::T => &'T',
            Self::N => &'N',
            Self::D => &'*',
        }
    }
}

impl Allele {
    pub fn is_known(&self) -> bool {
        match self {
            Self::N => false,
            _       => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn expected() -> HashMap<Allele, char> {
        HashMap::from_iter(vec![
            (Allele::A, 'A'),
            (Allele::C, 'C'),
            (Allele::G, 'G'),
            (Allele::T, 'T'),
            (Allele::N, 'N'),
        ])
    }
    #[test]
    fn display() {
        assert_eq!(format!("'{:_<5}'", Allele::A), "'A____'");
        assert_eq!(format!("'{:_<5}'", Allele::C), "'C____'");
        assert_eq!(format!("'{:_<5}'", Allele::G), "'G____'");
        assert_eq!(format!("'{:_<5}'", Allele::T), "'T____'");
        assert_eq!(format!("'{:_<5}'", Allele::N), "'N____'");
    }

    #[test]
    fn try_from_char() {
        for (allele, char) in expected() {
            assert_eq!(Allele::try_from(char).unwrap(), allele);
        }
    }

    #[test]
    #[should_panic]
    fn panic_try_from_char_panic() {
        Allele::try_from('x').expect("");
    }

}