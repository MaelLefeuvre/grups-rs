
use std::{fmt::{self, Display, Formatter}, ops::{Deref, DerefMut}, collections::HashMap};

use genome::Genome;

use itertools::Itertools;
use log::info;

use grups_io::write::GenericWriter;
use located_error::prelude::*;

pub mod pwd;
pub use pwd::Pwd;

mod individual;
pub use individual::Individual;

mod comparison;
pub use comparison::Comparison;

mod test;

mod error;
use error::ComparisonError;

const UNDEFINED_LABEL_PREFIX: &str = "Ind";
const DISPL_SEP             : &str  = " - ";
const PAIRS_FORMAT_LEN      : usize = 20;
const COUNT_FORMAT_LEN      : usize = 11;
const AVERG_FORMAT_LEN      : usize = 11;
const FLOAT_FORMAT_PRECISION: usize = 6;

//use super::{PAIRS_FORMAT_LEN, COUNT_FORMAT_LEN, AVERG_FORMAT_LEN, DISPL_SEP};


/// Vector of all the user-requested comparisons
#[derive(Debug, Default)]
pub struct Comparisons (Vec<Comparison>);


impl Display for Comparisons {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        self.0.iter().try_fold((), |_, comparison| {
            writeln!(f, "{comparison}")
        })
    }
}

impl Deref for Comparisons {
    type Target = [Comparison];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Comparisons {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


impl Comparisons {

    /// Instantiate and populate a new `Comparisons` object from the user-provided parameters.
    #[must_use]
    pub fn parse(
        individuals           : &[usize],
        min_depths            : &[u16],
        names                 : &[String],
        allow_self_comparison : bool,
        genome                : &Genome,
        blocksize             : u32
    ) -> Comparisons {
        let mut inds = vec![];
        for (i, index) in individuals.iter().enumerate() {
            let name = names.get(i);
            let min_depth = &min_depths[i % (min_depths.len())]; // wrap around min_depths if its length is lesser than the number of inds.
            inds.push(Individual::new(name, *index, *min_depth));
        }
        let mut comparisons: Vec<Comparison> = vec![];
        for pair in inds.iter().combinations_with_replacement(2) {
            let self_comparison = pair[0] == pair[1]; 
            if self_comparison && !allow_self_comparison {
                continue
            }
        comparisons.push(Comparison::new([pair[0].clone(), pair[1].clone()], self_comparison, genome, blocksize));            
        }
    Comparisons(comparisons)
    }

    
    /// Write the results of the observed average PWD in a `.pwd` output file.
    /// 
    /// # Errors
    /// - If the provided path to the output `.pwd` file is invalid or the user does not have the proper UNIX
    ///   permissions.
    /// - If any of the contained `Comparison` fails to get written to the output file. 
    pub fn write_pwd_results(&self, print_blocks: bool, output_files: &HashMap<String, String>) -> Result<()> {
        let mut pwd_writer = GenericWriter::new(Some(output_files["pwd"].clone()))?;
        
        let header = format!(
            "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$}{DISPL_SEP}\
             {: <AVERG_FORMAT_LEN$}",
             "Pair_name", "Raw.Overlap", "Raw.Sum.PWD", "Raw.Avg.PWD", "Raw.CI.95", "Raw.Avg.Phred");
        pwd_writer.write_iter(vec![&header])?; // Print PWD results to file.
        pwd_writer.write_iter(self.iter())?;   // 
        
        if log::log_enabled!(log::Level::Debug) {
            info!("Raw observed pairwise differences:");
            println!("{header}");
            println!("{self}"); // Print PWD results to console
        }

        if print_blocks {
            use genome::jackknife::{CHROM_FORMAT_LEN, COUNT_FORMAT_LEN, RANGE_FORMAT_LEN};
            let blk_header = format!(
            "{: <CHROM_FORMAT_LEN$}{DISPL_SEP}\
             {: <RANGE_FORMAT_LEN$}{DISPL_SEP}\
             {: <RANGE_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
             {: <COUNT_FORMAT_LEN$}",
             "chr", "start", "end", "overlap", "pwd"
            );
            for comparison in self.iter() {
                let pair = comparison.get_pair();
                let mut block_writer = GenericWriter::new(Some(output_files[&pair].clone()))?;
                block_writer.write_iter(vec![&blk_header])?;
                block_writer.write_iter(vec![&comparison.blocks])?;
            }
        }
        Ok(())
    }

    /// How many pairs are we comparing?
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Get the individual's pileup indices of each comparison
    #[must_use]
    pub fn get_pairs_indices(&self) -> Vec<[usize; 2]> {
        self.0.iter().map(Comparison::get_pair_indices).collect()
    }

    /// Check if this struct has been populated.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[must_use]
    pub fn get_pairs(&self) -> Vec<String> {
        self.0.iter().map(Comparison::get_pair).collect()
    }

    pub fn update_variance_unbiased(&mut self) {
        self.0.iter_mut().for_each(Comparison::update_variance_unbiased);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{fmt::Write, error::Error};
    use crate::comparisons::test::common;

    fn factorial(n: u32 ) -> u32 {
        match n {
            0   => 1,
            1.. => (1..n+1).product(),
        }
    }
    fn permutations_sample(n: u32) -> u32 {
        match n {
            0..=2 => factorial(n) - 1,
            3..   => factorial(n)/(factorial(n-2)*2),
        }
    }

    fn mock_comparisons(ind_set: &[usize], allow_self_comparison: bool) -> Comparisons {
        let min_depths = vec![2];
        let names = vec![];
        let genome = Genome::default();
        Comparisons::parse(ind_set, &min_depths, &names, allow_self_comparison, &genome, 50_000_000)
    }

    #[test]
    fn comparisons_length_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let comparisons = mock_comparisons(ind_set, true);
            let len = ind_set.len() as u32;
            assert_eq!(comparisons.len() as u32, permutations_sample(len)+len); 
        }
    }

    #[test]
    fn comparisons_length_no_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let comparisons = mock_comparisons(ind_set, false);
            let len = ind_set.len() as u32;
            assert_eq!(comparisons.len() as u32, permutations_sample(len)); 
        }
    }

    #[test]
    fn pair_indices_getter_self_allowed(){
        let ind_set = vec![0,1,2];
        let comparisons = mock_comparisons(&ind_set, true);
        let expected_output: Vec<Vec<usize>> = ind_set.into_iter().combinations_with_replacement(2).collect::<Vec<_>>();
        let results_vec: Vec<Vec<usize>> = comparisons.get_pairs_indices().iter().map(|x| x.to_vec()).collect();
        assert_eq!(results_vec, expected_output);
    }

    #[test]
    fn pair_indices_getter_no_self_allowed(){
        let ind_set = vec![0,1,2];
        let comparisons = mock_comparisons(&ind_set, false);
        let expected_output: Vec<Vec<usize>> = ind_set.into_iter().combinations(2).collect::<Vec<_>>();
        let results_vec: Vec<Vec<usize>> = comparisons.get_pairs_indices().iter().map(|x| x.to_vec()).collect();
        assert_eq!(results_vec, expected_output);
    }

    #[test]
    fn emptiness() {
        let mut empty_comp = Comparisons::default();
        assert!(empty_comp.is_empty());

        empty_comp.0.push(common::mock_comparison(false));
        assert!(! empty_comp.is_empty());
    }

    #[test]
    fn pair_string_getter_self_allowed(){
        let ind_set = vec![0,1,2];

        let mut expected_output = Vec::new();
        for indices in ind_set.iter().combinations_with_replacement(2) {
            let expected_label = format!("{UNDEFINED_LABEL_PREFIX}{}-{UNDEFINED_LABEL_PREFIX}{}", indices[0], indices[1]);
            expected_output.push(expected_label);
        }

        let comparisons = mock_comparisons(&ind_set, true);
        assert_eq!(comparisons.get_pairs(), expected_output);
    }

    #[test]
    fn pair_string_getter_no_self_allowed(){
        let ind_set = vec![0,1,2];

        let mut expected_output = Vec::new();
        for indices in ind_set.iter().combinations(2) {
            let expected_label = format!("{UNDEFINED_LABEL_PREFIX}{}-{UNDEFINED_LABEL_PREFIX}{}", indices[0], indices[1]);
            expected_output.push(expected_label);
        }

        let comparisons = mock_comparisons(&ind_set, false);
        assert_eq!(comparisons.get_pairs(), expected_output);
    }

    #[test]
    fn display() -> Result<(), Box<dyn Error>> {
        let ind_set = vec![0,1,2];
        let comparisons = mock_comparisons(&ind_set, false);
        let mut expected_output = String::new();
        for indices in ind_set.iter().combinations(2) {
            let expected_pair_name=format!("{UNDEFINED_LABEL_PREFIX}{}-{UNDEFINED_LABEL_PREFIX}{}", indices[0], indices[1]);
            writeln!(expected_output,
                "{: <PAIRS_FORMAT_LEN$}{DISPL_SEP}\
                {: <COUNT_FORMAT_LEN$}{DISPL_SEP}\
                {: <AVERG_FORMAT_LEN$.1}{DISPL_SEP}\
                {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
                {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
                {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
                expected_pair_name, 0, 0.0, f64::NAN, f64::NAN, f64::NAN
            )?;
        }

        assert_eq!(format!("{comparisons}"), expected_output);
        Ok(())
    }

}