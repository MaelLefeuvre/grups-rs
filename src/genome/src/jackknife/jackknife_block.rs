
use crate::snpcoord::SNPCoord;
use std::hash::{Hash, Hasher};
use std::ops::Range;
use std::fmt;

use super::{CHROM_FORMAT_LEN, COUNT_FORMAT_LEN, RANGE_FORMAT_LEN};

pub struct Pseudovalue {
    pub hj     : f64,
    pub theta_j: f64,
}

impl Pseudovalue {
    pub fn weigthed_pseudovalue(&self) -> f64 {
        self.theta_j / self.hj
    }
}

///Chromosome Block for Jackknife resampling. Implemented within struct `JackknifeBlocks`, which is itself implemented
///within structs `pwd_from_stdin::pileup::Comparison`
/// - chromosome  : chromosome name (warn: not by index -> cannot handle specific names, such as "X", "MT", "chr7", etc.)
/// - range       : [start, end[ coordinates of the block.
/// - site_counts : counts the number of overlapping SNPs for a given pair of individuals.
/// - pwd_counts  : counts the number of pairwise differences for a given pair of individuals.
/// 
/// # Traits :
///  - `Eq` and `PartialEq` : with `Self`.
///  - `Eq` and `PartialEq` : with struct `pwd_from_stdin::genome::SNPCoord`
///     - mainly implemented in order to retrieve the block corresponding with a given SNP position.
///       See: `Self::find_block()` or `JackknifeBlocks::find_blocks()`
///  - `Hash`               : along `chromosome` and `range` fields.
///  - `Display`            : Pretty print for file and/or console output Recursively called by `JackknifeBlocks` when it itself is displayed.
#[derive(Debug)]
pub struct JackknifeBlock {
    pub chromosome : u8,
    pub range      : Range<u32>,
    site_counts: u32,
    pwd_counts : u32,
}

impl JackknifeBlock {
    pub fn new(chromosome: u8, start: u32, end: u32) -> JackknifeBlock {
        JackknifeBlock{chromosome, range: Range{start, end}, site_counts: 0, pwd_counts:0}
    }

    /// Incrementer for the `pwd_counts` field. Called by `pwd_from_stdin::Comparison::compare()`
    /// when a pairwise difference has been declared.
    pub fn add_pwd(&mut self) {
        self.pwd_counts += 1;
    }

    /// Incremented for the `site_counts` field. Infaillibly called by `pwd_from_stdin::Comparison::compare()`
    pub fn add_count(&mut self) {
        self.site_counts += 1;
    }

    pub fn compute_unequal_delete_m_pseudo_value(&self, sum_pwd: u32, sum_overlap: u32) -> Pseudovalue {
        // hj = n / m_j
        // theta: Estimate of avg_PWD bassed on all the observations.
        // theta_minus_j: Estimate of avg_PWD based on all the observations, except those from Block j
        let hj          = sum_overlap as f64 / self.site_counts as f64; // hj = n/m_j
        let theta         = sum_pwd as f64 / sum_overlap as f64;          // Estimate of avg_pwd based on all the observations 
        let theta_minus_j = (sum_pwd - self.pwd_counts) as f64 / (sum_overlap - self.site_counts) as f64; // Estimate of avg_pwd based on all the observations 

        let theta_j = hj * theta - (hj - 1.0) * theta_minus_j;
        Pseudovalue{hj, theta_j}
    }
}

impl fmt::Display for JackknifeBlock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
            "{: <CHROM_FORMAT_LEN$} - \
             {: <RANGE_FORMAT_LEN$} - \
             {: <RANGE_FORMAT_LEN$} - \
             {: <COUNT_FORMAT_LEN$} - \
             {: <COUNT_FORMAT_LEN$}",
            self.chromosome,
            self.range.start,
            self.range.end,
            self.site_counts,
            self.pwd_counts
        )
    }
}

impl PartialEq<SNPCoord> for JackknifeBlock {
    fn eq(&self, other: &SNPCoord) -> bool {
        self.chromosome == other.chromosome && self.range.start <= other.position && self.range.end > other.position
    }
}

impl PartialEq<JackknifeBlock> for JackknifeBlock {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.range == other.range 
    }
}

impl Eq for JackknifeBlock {}

impl Hash for JackknifeBlock {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.range.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const N_ITERS: u32 = 1_000_000;

    #[test]
    fn increment_pwd_count() {
        let mut block = JackknifeBlock::new(1, 0, 1000);
        assert_eq!(block.pwd_counts, 0);
        block.add_pwd();
        assert_eq!(block.pwd_counts, 1);
    }

    #[test]
    fn increment_site_count() {
        let mut block = JackknifeBlock::new(1, 0, 1000);
        assert_eq!(block.site_counts, 0);
        block.add_count();
        assert_eq!(block.site_counts, 1);
    }

    #[test]
    fn display(){
        let block = JackknifeBlock::new(1, 0, 1000);
        let expect_out = format!(
            "{: <CHROM_FORMAT_LEN$} - \
             {: <RANGE_FORMAT_LEN$} - \
             {: <RANGE_FORMAT_LEN$} - \
             {: <COUNT_FORMAT_LEN$} - \
             {: <COUNT_FORMAT_LEN$}",
            1, 0, 1000, 0, 0
        );
        assert_eq!(expect_out, format!("{block}"));
    }

    #[test]
    fn snpcoord_equality() {
        let block = JackknifeBlock::new(10, 11_000, 12_000);

        //snpcoord == block if block.start <= snpcoord.position < block.end AND block.chr == snpcoord.chr
        assert_eq!(block, SNPCoord{chromosome: 10, position: 11_000, reference: None, alternate: None}); // le. range    ; same chromosome.
        assert_eq!(block, SNPCoord{chromosome: 10, position: 11_500, reference: None, alternate: None}); // within range ; same chromosome.
    }

    #[test]
    fn snpcoord_inequality() {
        let block = JackknifeBlock::new(10, 11_000, 12_000);

        assert_ne!(block, SNPCoord{chromosome: 10, position: 10_999, reference: None, alternate: None}); // lt. range    ; same chromosome
        assert_ne!(block, SNPCoord{chromosome: 10, position: 12_000, reference: None, alternate: None}); // ge. range    ; same chromosome
        assert_ne!(block, SNPCoord{chromosome: 10, position: 12_001, reference: None, alternate: None}); // gt. range    ; same chromosome
        assert_ne!(block, SNPCoord{chromosome: 11, position: 11_500, reference: None, alternate: None}); // within range ; different chromosome.
    }

    #[test]
    fn block_equality() {
        //        self.chromosome == other.chromosome && self.range == other.range 
        let block_0 = JackknifeBlock::new(10, 11_000, 12_000);
        let block_1 = JackknifeBlock::new(10, 11_000, 12_000);
        assert_eq!(block_0, block_1)
    }

    #[test]
    fn block_inequality() {
        //        self.chromosome == other.chromosome && self.range == other.range 
        let (chr, start, end) = (10, 11_000, 12_000);
        let block = JackknifeBlock::new(10, 11_000, 12_000);

        let deviations: [i32; 3] = [-1, 0, 1];
        for chr_deviation in deviations {
            for start_deviation in deviations {
                for end_deviation in deviations {
                    let chr   = chr   + chr_deviation;
                    let start = start + start_deviation;
                    let end   = end   + end_deviation;

                    let other_block = JackknifeBlock::new(chr as u8, start as u32, end as u32);
                    if (chr_deviation,start_deviation,end_deviation) == (0,0,0) { // Only case where we expect equality
                        assert_eq!(block, other_block);
                    } else {
                        assert_ne!(block, other_block);
                    }
                }
            }
        }
    }

    #[test]
    fn hash_block() {
        let mut test_hashset = std::collections::HashSet::new();
        let ranges: Vec<u32> = (1000..N_ITERS).step_by(1000).collect();
        let chr = 10;

        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);
            let block = JackknifeBlock::new(chr, start, end);
            assert!(test_hashset.insert(block));
        }

        // Retrieve same blocks
        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);

            let deviations: [i32; 3] = [-1, 0, 1];
            for chr_deviation in deviations {
                for start_deviation in deviations {
                    for end_deviation in deviations {
                        let chr   = chr   as i32 + chr_deviation;
                        let start = start as i32 + start_deviation;
                        let end   = end   as i32 + end_deviation;
                        let other_block = JackknifeBlock::new(chr as u8, start as u32, end as u32);

                        if (chr_deviation,start_deviation,end_deviation) == (0,0,0) { // Only  case where we expect equality
                            assert!( test_hashset.contains(&other_block) );
                        } else {
                            assert!( !test_hashset.contains(&other_block) );
                        }
                    }
                }
            }
        }

    }

}