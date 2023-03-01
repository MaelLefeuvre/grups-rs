use std::{fmt, collections::BTreeSet};

use genome::jackknife::{JackknifeBlocks, JackknifeEstimates};
use genome::Genome;
use located_error::LocatedError;

use crate::pileup::{Pileup, Line};
use super::ComparisonError;
use super::{Individual, Pwd};
use super::{PAIRS_FORMAT_LEN, COUNT_FORMAT_LEN, AVERG_FORMAT_LEN, DISPL_SEP, FLOAT_FORMAT_PRECISION};

use anyhow::Result;


use log::warn;

// One pass standard deviation calculator: http://suave_skola.varak.net/proj2/stddev.pdf
#[derive(Debug)]
struct Variance{
    meansum: f64,
    std_dev_sum: f64,
    n: usize,
}

impl Variance {
    pub fn new() -> Self {
        Self{meansum: 0.0, std_dev_sum:0.0, n:0}
    }

    /// Unbiased, two-pass variance estimation algorithm 
    pub fn update_from_iter(&mut self, avg_pwd: f64, pwds: &BTreeSet<Pwd>) {
        self.std_dev_sum = pwds.iter()
            .map(|pwd| (avg_pwd - pwd.pwd).powf(2.0))
            .sum::<f64>();

        self.n = pwds.len();

    }

    /// Single-pass variance estimation algorithm.
    #[allow(dead_code)]
    pub fn dynamic_update(&mut self, pwd: f64) {
        if self.n == 0 {
            self.meansum = pwd;
        };
        self.n += 1;
        let stepsum = pwd - self.meansum;
        let stepmean = ((self.n - 1) as f64 * stepsum) / self.n as f64;
        self.meansum += stepmean;
        self.std_dev_sum += stepmean * stepsum;
    }

    pub fn std_dev(&self) -> f64 {
        if self. n == 0 {
            return 0.0
        }
        (self.std_dev_sum / (self.n as f64 - 1.0)).sqrt()
    }

    pub fn confidence_interval(&self) -> f64 {
        1.96 * (self.std_dev()/(self.n as f64).sqrt())
    }
}

/// A struct representing a given pairwise estimation of relatedness between two individuals.
/// - `pair`            : contains a representation of the individuals being compared.
/// - `self_comparison` : whether or not this comparison is a Self-comparison. (i.e. pair.0 == pair.1)
/// - `pwd`             : counter for the number of pairwise differences found between our pair.
/// - `sum_phred`       : sum of the avg. phred-scores of each overlapping SNP. mainly used to compute the average
///                       Phred-score (i.e. `overlap`/`sum_phred`)
/// - `blocks`          : genome blocks used for jackknife resampling.
/// 
/// # Traits : `Debug`
/// 
/// # @TODO:
///  - `variance` and `positions` should be tied within their own struct.
#[derive(Debug)]
pub struct Comparison {
    pair            : [Individual; 2],
    self_comparison : bool,
    variance        : Variance,
    pub blocks      : JackknifeBlocks,
    pub positions   : BTreeSet<Pwd>
}

impl Comparison {
    pub fn new(mut pair: [Individual; 2], self_comparison: bool, genome: &Genome, blocksize: u32) -> Comparison {
        if self_comparison {
            let pair_name = format!("{}-{}", pair[0].name, pair[1].name);
            for individual in &mut pair {
                if individual.min_depth < 2 {
                    warn!("[{pair_name}]: A minimal depth of 2 is required for self comparisons, while {} has min_depth = {}.\n
                        Rescaling said value to 2 (for this specific comparison)...", individual.name, individual.min_depth
                    );
                    individual.min_depth = 2;
                }
            }
        }
        Comparison {pair, self_comparison, variance: Variance::new(), blocks: JackknifeBlocks::new(genome, blocksize), positions: BTreeSet::new()}
    }

    pub fn get_pair_indices(&self) -> [usize; 2] {
        [self.pair[0].index, self.pair[1].index]
    }

    /// Check if the sequencing depth of a given overlap is over the minimum required sequencing depth for each individual.
    pub fn satisfiable_depth(&self, pileups: &[Pileup]) -> bool {
        self.pair[0].satisfiable_depth(pileups) && self.pair[1].satisfiable_depth(pileups)
    }

    /// Compare our two individuals at the given SNP position ; increment the appropriate counters after the comparison 
    /// has been made.
    pub fn compare(&mut self, line: &Line) -> Result<()> {
        //let pwd = Pwd::one(line.coordinate, &random_nucl);
        let loc_msg = || format!("While comparing pair {:?}", self.pair);
        let pwd = match self.self_comparison {
            true  => Pwd::deterministic_self(line, &self.pair),
            false => Pwd::deterministic_pairwise(line, &self.pair)
        };

        let current_block = self.blocks
            .find_block(&line.coordinate)
            .ok_or(ComparisonError::MissingBlock(line.coordinate)).with_loc(loc_msg)?;

        current_block.add_count();
        current_block.add_pwd(pwd.avg_local_pwd());

        self.positions.insert(pwd);
        Ok(())
    }

    pub fn get_sum_pwd(&self) -> f64 {
        self.positions.iter()
            .map(Pwd::avg_local_pwd)
            .sum::<f64>()
    }

    // Getter for the average pairwise difference across our overlapping snps.
    pub fn get_avg_pwd(&self) -> f64 {
        self.get_sum_pwd() / self.positions.len() as f64
    }

    pub fn get_overlap(&self) -> usize {
        self.positions.len()
    }

    fn get_sum_phred(&self) -> f64 {
        self.positions.iter()
            .map(Pwd::compute_avg_phred)
            .sum::<f64>()
    }
    // Getter for the average pairwise phred score across our overlapping snps.
    pub fn get_avg_phred(&self) -> f64 {
        self.get_sum_phred() / self.positions.len() as f64
    }

    // Return a formated string representing our pair of individuals. 
    pub fn get_pair (&self,)-> String {
        format!("{}-{}", &self.pair[0].name, &self.pair[1].name)
    }

    /// Obtain the jacknife estimates of the avg. PWD's standard deviation at each block.
    pub fn get_jackknife_estimates(&self) -> JackknifeEstimates {
        self.blocks.compute_unequal_delete_m_pseudo_values(self.get_sum_pwd(), self.positions.len() as u32)
    }

    /// Optain the 95% confidence interval for the observed avg. PWD
    pub fn get_confidence_interval(&self) -> f64 {
        self.variance.confidence_interval()
    }

    /// Compute the variance and standard deviation for the observed avg. PWD using a two-pass method.
    pub fn update_variance_unbiased(&mut self) {
        self.variance.update_from_iter(self.get_avg_pwd(), &self.positions);
    }

    //// WIP: heterozygocity ratio
    //pub fn get_heterozygocity_ratio(&self) -> f64 {
    //    let sum_hom_alt = self.positions.iter().map(|pwd| pwd.hom_alt_sum).sum::<f64>();
    //    self.get_sum_pwd() / sum_hom_alt
    //}

}

impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
            "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.1}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
            self.get_pair(),
            self.positions.len(),
            self.get_sum_pwd(),
            self.get_avg_pwd(),
            self.get_confidence_interval(),
            self.get_avg_phred()
        )
    }
}

#[cfg(test)]
mod tests {
    use genome::SNPCoord;

    use crate::pileup;
    use crate::comparisons::test::common;
    use super::*;
    use super::super::UNDEFINED_LABEL_PREFIX;


    #[test]
    fn pair_indices_getter() {
        let comparison = common::mock_comparison(false);
        assert_eq!(comparison.get_pair_indices(), common::MOCK_IND_INDICES)
    }

    #[test]
    fn satisfiable_depth_both() -> Result<()> {
        // If both individual has depth >= pileup.depth ==> return false. 
        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[2,2], 30, false)?;
        assert!(comparison.satisfiable_depth(&pileups[..]));
        Ok(())
    }

    #[test]
    fn satisfiable_depth_one() -> Result<()> {
        // If one or the other individual has depth < pileup.depth ==> return false. 
        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[3,1], 30, false)?;
        assert!( !comparison.satisfiable_depth(&pileups[..]));

        let pileups = common::mock_pileups(&[1,3], 30, false)?;
        assert!( !comparison.satisfiable_depth(&pileups[..]));
        Ok(())
    }

    #[test]
    fn satisfiable_depth_none() -> Result<()> {
        // If both individual has depth < pileup.depth ==> return false. 
        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[1,1], 30, false)?;
        assert!( !comparison.satisfiable_depth(&pileups[..]));
        Ok(())
    }


    fn test_compare(
        comparison: &mut Comparison,
        positions : &[u32],
        ref_base  : char,
        nucs      : [&str; 2],
        quals     : [&str; 2]
    ) -> Result<()> {
        for pos in positions {
            let raw_line=format!("22\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                *pos,
                ref_base,
                nucs[0].len(), nucs[0], quals[0],
                nucs[1].len(), nucs[1], quals[1])
            ;
            let line = pileup::Line::new(raw_line.as_str(), true)?;
            comparison.compare(&line)?;
        }
        Ok(())
    }

    #[test]
    fn test_compare_pairwise_no_pwd() -> Result<()> {
        let mut mock_comparison = common::mock_comparison(false);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "TT"], ["JJ", "AA"])?;

        assert_eq!(mock_comparison.get_sum_phred(), 73.0);       // (J  +  A)
        assert_eq!(mock_comparison.get_avg_phred(), 73.0 / 2.0);
        assert_eq!(mock_comparison.get_sum_pwd(), 0.0);          // (T <-> T) * 2
        assert_eq!(mock_comparison.get_avg_pwd(), 0.0);
        assert_eq!(mock_comparison.positions.len(), 2);
        Ok(())
    }

    #[test]
    fn test_compare_self_no_pwd() -> Result<()> {
        let mut mock_comparison = common::mock_comparison(true);

        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "GG"], ["JJ", "AA"])?;

        assert_eq!(mock_comparison.get_sum_phred(), 82.0);        // (J + *J*) (NOT A, since we're self-comparing)
        assert_eq!(mock_comparison.get_avg_phred(), 82.0 / 2.0);
        assert_eq!(mock_comparison.get_sum_pwd(), 0.0);           // (T <-> *T*) * 2 (NOT G, since we're self-comparing)
        assert_eq!(mock_comparison.get_avg_pwd(), 0.0);          
        assert_eq!(mock_comparison.positions.len(), 2);         
        Ok(())
    }


    #[test]
    fn test_compare_pairwise_pwd() -> Result<()> {
        let mut mock_comparison = common::mock_comparison(false);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "CC"], ["JJ", "AA"])?;

        assert_eq!(mock_comparison.get_sum_phred(), 73.0);       // (J + A)
        assert_eq!(mock_comparison.get_avg_phred(), 73.0 / 2.0);
        assert_eq!(mock_comparison.get_sum_pwd(), 2.0);          // (T <-> C) *2
        assert_eq!(mock_comparison.get_avg_pwd(), 1.0);
        assert_eq!(mock_comparison.positions.len(), 2);
        Ok(())
    }

    #[test]
    fn test_phred_increment() -> Result<()> {
        let mut mock_comparison = common::mock_comparison(false);
        let mut expected_sum_phred: u32 = 0;
        let test_quals = ["5", "?"];
        let avg_test_qual = 25; // '5' + '?' /2
        for pos in 1..1000u32 {
            test_compare(&mut mock_comparison, &[pos], 'C', ["T", "T"], test_quals)?;


            // comparison.sum_phred represent the average of two given nucleotide's scores. i.e. (J + A) / 2
            expected_sum_phred += avg_test_qual; 

            // Ensure comparison.sum_phred represent the sum average of two given nucleotide's scores
            assert_eq!(mock_comparison.get_sum_phred(), expected_sum_phred as f64);

            assert_eq!(mock_comparison.get_avg_phred(), avg_test_qual as f64); // Avg. pwd should stay the same. 
        }
        Ok(())
    }

    #[test]
    fn update_positions() -> Result<()> {
        let mut mock_comparison = common::mock_comparison(false);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "TT"], ["JJ", "JJ"])?;
        test_compare(&mut mock_comparison, &[30,40], 'C', ["TT", "CC"], ["JJ", "AA"])?; 
        // J: 41 | A: 32
        assert_eq!(mock_comparison.get_sum_phred(), 155.0);       // (J + A) 
        assert_eq!(mock_comparison.get_avg_phred(), 155.0 / 4.0);
        assert_eq!(mock_comparison.get_sum_pwd(), 2.0);           // (T <-> C) *2
        assert_eq!(mock_comparison.get_avg_pwd(), 0.5);
        assert_eq!(mock_comparison.positions.len(), 4);

        
        // Remove the last position from the comparisons.
        let x = SNPCoord::try_new(22, 40, 'N', 'N')?;
        mock_comparison.positions.remove(&Pwd::initialize(x.coordinate));

        // Check pwd statistics are updated and sane...
        assert_eq!(mock_comparison.get_sum_phred(), 118.5);       // (J + A)
        assert_eq!(mock_comparison.get_avg_phred(), 118.5 / 3.0);
        assert_eq!(mock_comparison.get_sum_pwd(), 1.0);           // (T <-> C) *2
        assert_eq!(mock_comparison.get_avg_pwd(), 0.3333333333333333);
        assert_eq!(mock_comparison.positions.len(), 3);

        Ok(())
    }


    #[test]
    fn display() {
        let mock_comparison = common::mock_comparison(false);
        let expected_pair_name=format!("{UNDEFINED_LABEL_PREFIX}0-{UNDEFINED_LABEL_PREFIX}1");
        let expect_out = format!(
            "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.1}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
             expected_pair_name, 0, 0.0, f64::NAN, f64::NAN, f64::NAN
        );
        assert_eq!(expect_out, format!("{mock_comparison}"));
    }
}
