use std::{cmp::Ordering, str::FromStr};

use log::debug;
use located_error::LocatedError;
use crate::coordinate::ChrIdx;

mod error;
use error::ChromosomeError;

/// A simple struct representing a chromosome. This is mainly used to compute Jackknife Blocks
/// # Fields:
/// - `index` : 0-based index of the chromosome
/// - `name`  : 1-based name of the chromosome (i.e. 1-22)
/// - `length`: length of the chromosome (bp)
#[derive(Debug, Clone, Copy)]
pub struct Chromosome {
    //index  : usize,
    pub name   : ChrIdx,
    pub length : u32,
}

impl Chromosome {
    /// Instantiate a new chromosome
    /// # Fields:
    /// - `index` : 0-based index of the chromosome
    /// - `name`  : 1-based name of the chromosome (i.e. 1-22)
    /// - `length`: length of the chromosome (bp)
    #[must_use]
    pub fn new(name: impl Into<ChrIdx>, length: u32) -> Chromosome{
        Chromosome{name: name.into(), length}
    }
}

impl PartialEq for Chromosome {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}
impl Eq for Chromosome {}

impl PartialOrd for Chromosome {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chromosome {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.name).cmp(&other.name)
    }
}


impl FromStr for Chromosome {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use ChromosomeError::*;
        let split_line: Vec<&str> = s.split('\t').collect();
        let name = ChrIdx::from_str(split_line[0])
            .with_loc(|| ParseChrIdx)?;
        let length = split_line[1].parse::<u32>()
            .with_loc(||ParseLength)?;

        debug!("Chromosome: {: <10} {: <12}", name, length);
        Ok(Self::new(name, length))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_two_chromosome(name_first: u8, name_second: u8) -> [Chromosome; 2] {
        let chr1 = Chromosome::new(name_first, 1000);
        let chr2 = Chromosome::new(name_second, 2000);
        [chr1, chr2]
    }

    #[test]
    fn chromosome_name_equality() {
        let test_chr = get_two_chromosome(1, 1);
        assert_eq!(test_chr[0], test_chr[1]);
    }

    #[test]
    fn chromosome_name_inequality() {
        let test_chr = get_two_chromosome(1, 22);
        assert!(test_chr[0] != test_chr[1]);
    }

    #[test]
    fn chromosome_name_lge(){
        let test_chr = get_two_chromosome(1, 1);
        assert!(test_chr[0] <= test_chr[1]);
        assert!(test_chr[0] >= test_chr[1]);

    }

    #[test]
    fn chromosome_name_lgt(){
        let test_chr = get_two_chromosome(1, 2);
        assert!(test_chr[0] < test_chr[1]);
        assert!(test_chr[1] > test_chr[0]);
    }

    #[test]
    fn from_str_ok() {
        for i in 1..=22 {
            let chr = Chromosome::from_str(&format!("{i}\t555555"));
            assert!(chr.is_ok())
        }
    }

    #[test]
    fn from_str_missing_chr(){
        let chr = Chromosome::from_str("55555");
        assert!(chr.is_err());
    }

    #[test]
    fn from_str_missing_length(){
        let chr = Chromosome::from_str("22\t");
        assert!(chr.is_err());
    }

    #[test]
    fn from_str_xchr() {
        for label in ["X", "chrX"].iter() {
            let chr = Chromosome::from_str(&format!("{label}\t123456789"));
            assert!(chr.is_ok_and(|c| c.name == ChrIdx(88)));
        }
    }
}