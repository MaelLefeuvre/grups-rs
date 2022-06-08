use genome::jackknife::{JackknifeBlocks, JackknifeEstimates};
use genome::{SNPCoord, Genome};
use crate::pileup::{Pileup, Line, Nucleotide};
use std::error::Error;
use std::{
    fmt,
    collections::BTreeSet
};

use super::{Individual, Pwd};
use super::{PAIRS_FORMAT_LEN, COUNT_FORMAT_LEN, AVERG_FORMAT_LEN, DISPL_SEP, FLOAT_FORMAT_PRECISION};
/// A struct representing a given pairwise estimation of relatedness between two individuals.
/// - pair            : contains a representation of the individuals being compared.
/// - self_comparison : whether or not this comparison is a self_comparison. (i.e. pair.0 == pair.1)
/// - overlap         : counter for the number of overlapping SNP positions between our pair.
/// - pwd             : counter for the number of pairwise differences found between our pair.
/// - sum_phred       : sum of the avg. phred-scores of each overlapping SNP. mainly used to compute the average
///                     Phred-score (i.e. overlap/sum_phred)
/// - blocks          : genome blocks used for jackknife resampling.
/// 
/// # Traits : `Debug`
#[derive(Debug)]
pub struct Comparison {
    pair            : [Individual; 2],
    self_comparison : bool,
    pwd             : u32,
    sum_phred       : u32,
    pub blocks      : JackknifeBlocks,
    pub positions   : BTreeSet<Pwd>
}

impl Comparison {
    pub fn new(pair: [Individual; 2], self_comparison: bool, genome: &Genome, blocksize: u32) -> Comparison {
        Comparison {pair, self_comparison, pwd:0, sum_phred:0, blocks: JackknifeBlocks::new(genome, blocksize), positions: BTreeSet::new()}
    }

    pub fn get_pair_indices(&self) -> [usize; 2] {
        [self.pair[0].index, self.pair[1].index]
    }

    // Check if the sequencing depth of a given overlap is over the minimum required sequencing depth for each individual.
    pub fn satisfiable_depth(&self, pileups: &[Pileup]) -> bool {
        self.pair[0].satisfiable_depth(pileups) && self.pair[1].satisfiable_depth(pileups)
    }

    // Compare our two individuals at the given SNPposition ; increment the appropriate counters after the comparison 
    // has been made.
    pub fn compare(&mut self, line: &Line) -> Result<(), Box<dyn Error>> {
        //self.positions.insert(line.coordinate.clone());

        let random_nucl: Vec<&Nucleotide> = if self.self_comparison {   // Self comparison => Combination without replacement. 
            line.random_sample_self(&self.pair[0].index)                //   - possible combinations: n!/(n-2)!
        } else {                                                        // Std comparison  => Permutation.
            line.random_sample_pair(&self.pair)                         //  - possible combinations: n_1 * n_2
        };

        let pwd = Pwd::new(line.coordinate, &random_nucl);
        
        self.add_phred(&random_nucl);

        let current_block = self.blocks
            .find_block(&line.coordinate)
            .ok_or(format!("Cannot find corresponding Jackknife Block for {}", line.coordinate))?;

        current_block.add_count();
        if pwd.is_pwd() {
            self.pwd +=1;
            current_block.add_pwd();
        }

        self.positions.insert(pwd);
        Ok(())
    }

    // Increment our `sum_phred` counter. The incremented value is computed as the average of the two sampled nucleotides.
    fn add_phred<'a>(&mut self, nuc: &[&'a Nucleotide]) {
        self.sum_phred += ( (nuc[0].phred + nuc[1].phred)/2 ) as u32;
    }

    // Check if there is a pairwise difference.
    fn check_pwd<'a>(nuc: &[&'a Nucleotide]) -> bool {
        nuc[0].base != nuc[1].base
    }

    // Getter for the average pairwise difference across our overlapping snps.
    pub fn get_avg_pwd(&self) -> f64 {
        self.pwd as f64 / self.positions.len() as f64
    }

    // Getter for the average pairwise phred score across our overlapping snps.
    pub fn get_avg_phred(&self) -> f64 {
        self.sum_phred as f64/self.positions.len() as f64
    }

    // Return a formated string representing our pair of individuals. 
    pub fn get_pair (&self,)-> String {
        format!("{}-{}", &self.pair[0].name, &self.pair[1].name)
    }

    pub fn get_jackknife_estimates(&self) -> JackknifeEstimates {
        self.blocks.compute_unequal_delete_m_pseudo_values(self.pwd, self.positions.len() as u32)
    }
}


impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let jackknife_estimates = self.get_jackknife_estimates();
        write!(f,
            "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
            self.get_pair(),
            self.positions.len(),
            self.pwd,
            self.get_avg_pwd(),
            jackknife_estimates.variance,
            self.get_avg_phred()
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::pileup;
    use crate::comparisons::tests::common;
    use super::*;
    use super::super::UNDEFINED_LABEL_PREFIX;

    #[test]
    fn pair_indices_getter() {
        let comparison = common::mock_comparison(false);
        assert_eq!(comparison.get_pair_indices(), common::MOCK_IND_INDICES)
    }

    #[test]
    fn satisfiable_depth_both() {
        // If both individual has depth >= pileup.depth ==> return false. 

        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[2,2], 30, false);
        assert!(comparison.satisfiable_depth(&pileups[..]));
    }

    #[test]
    fn satisfiable_depth_one() {
        // If one or the other individual has depth < pileup.depth ==> return false. 
        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[3,1], 30, false);
        assert!( !comparison.satisfiable_depth(&pileups[..]));

        let pileups = common::mock_pileups(&[1,3], 30, false);
        assert!( !comparison.satisfiable_depth(&pileups[..]));

    }

    #[test]
    fn satisfiable_depth_none() {
        // If both individual has depth < pileup.depth ==> return false. 
        let comparison = common::mock_comparison(false);
        let pileups = common::mock_pileups(&[1,1], 30, false);
        assert!( !comparison.satisfiable_depth(&pileups[..]));
    }


    fn test_compare(comparison: &mut Comparison, positions: &[u32], ref_base: char, nucs: [&str; 2], quals: [&str; 2]) {
        for pos in positions {
            let raw_line=format!("22\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                *pos,
                ref_base,
                nucs[0].len(), nucs[0], quals[0],
                nucs[1].len(), nucs[1], quals[1])
            ;
            let line = pileup::Line::new(&raw_line.as_str(), true).unwrap();
            comparison.compare(&line).unwrap();
        }
    }

    #[test]
    fn test_compare_pairwise_no_pwd() {
        let mut mock_comparison = common::mock_comparison(false);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "TT"], ["JJ", "AA"]);

        assert_eq!(mock_comparison.sum_phred, 72);                 // ('J'  +  'A')
        assert_eq!(mock_comparison.get_avg_phred(), 72.0 / 2.0);
        assert_eq!(mock_comparison.pwd, 0);                        // ('T' <-> 'T') * 2
        assert_eq!(mock_comparison.get_avg_pwd(), 0.0);
        assert_eq!(mock_comparison.positions.len(), 2);
    }

    #[test]
    fn test_compare_self_no_pwd() {
        let mut mock_comparison = common::mock_comparison(true);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "TT"], ["JJ", "AA"]);

        assert_eq!(mock_comparison.sum_phred, 82);               // ('J'  +  'J')
        assert_eq!(mock_comparison.get_avg_phred(), 82.0 / 2.0);
        assert_eq!(mock_comparison.pwd, 0);                      // ('T' <-> 'T') * 2
        assert_eq!(mock_comparison.get_avg_pwd(), 0.0);          
        assert_eq!(mock_comparison.positions.len(), 2);         
    }


    #[test]
    fn test_compare_pairwise_pwd() {
        let mut mock_comparison = common::mock_comparison(false);
        test_compare(&mut mock_comparison, &[10,20], 'C', ["TT", "CC"], ["JJ", "AA"]);

        assert_eq!(mock_comparison.sum_phred, 72);               // ('J'  +  'A')
        assert_eq!(mock_comparison.get_avg_phred(), 72.0 / 2.0);
        assert_eq!(mock_comparison.pwd, 2);                      // ('T' <-> 'C')
        assert_eq!(mock_comparison.get_avg_pwd(), 1.0);
        assert_eq!(mock_comparison.positions.len(), 2);
    }

    #[test]
    fn test_phred_increment() {
        let mut mock_comparison = common::mock_comparison(false);
        let mut expected_sum_phred: u32 = 0;
        let test_quals = ["5", "?"];
        let avg_test_qual = 25; // '5' + '?' /2
        for pos in 1..1000u32 {
            test_compare(&mut mock_comparison, &[pos], 'C', ["T", "T"], test_quals);


            // comparison.sum_phred represent the average of two given nucleotide's scores. i.e. (J + A) / 2
            expected_sum_phred += avg_test_qual; 

            // Ensure comparison.sum_phred represent the sum average of two given nucleotide's scores
            assert_eq!(mock_comparison.sum_phred, expected_sum_phred);

            assert_eq!(mock_comparison.get_avg_phred(), avg_test_qual as f64); // Avg. pwd should stay the same. 
        }
    }

    #[test]
    fn pwd_checker_identity() {
        let nucleotides = [
            &Nucleotide{base:'A', phred: 20},
            &Nucleotide{base:'A', phred: 30}
        ];
        assert_eq!(Comparison::check_pwd(&nucleotides[..]), false);
    }

    #[test]
    fn pwd_checker_mismatch() {
        let nucleotides = [
            &Nucleotide{base:'C', phred: 20},
            &Nucleotide{base:'A', phred: 30}
        ];
        assert_eq!(Comparison::check_pwd(&nucleotides[..]), true);
    }

    #[test]
    fn display() {
        let mock_comparison = common::mock_comparison(false);
        let expected_pair_name=format!("{UNDEFINED_LABEL_PREFIX}0-{UNDEFINED_LABEL_PREFIX}1");
        let expect_out = format!(
            "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
             expected_pair_name, 0, 0, f64::NAN, f64::NAN
        );
        assert_eq!(expect_out, format!("{mock_comparison}"));
    }
}
