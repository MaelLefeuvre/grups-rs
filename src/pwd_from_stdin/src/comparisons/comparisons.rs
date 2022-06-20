use genome::Genome;
use std::fmt;
use std::ops::{Deref, DerefMut};

use itertools::Itertools;

use super::comparison::Comparison;
use super::individual::Individual;

/// Vector of all the user-requested comparisons
#[derive(Debug, Default)]
pub struct Comparisons (Vec<Comparison>);

impl Comparisons {
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
            let min_depth = &min_depths[(i % (min_depths.len())) as usize]; // wrap around min_depths if its length is lesser than the number of inds.
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

    ///How many pairs are we comparing? 
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn get_pairs_indices(&self) -> Vec<[usize; 2]> {
        Vec::from_iter(self.0.iter().map(|comparison| comparison.get_pair_indices()))
    }

    // Clippy is complaining when i'm not implementing this public method :)
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn get_pairs(&self) -> Vec<String> {
        self.0.iter().map(|c| c.get_pair()).collect()
    }
}

impl fmt::Display for Comparisons {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        self.0.iter().fold(Ok(()), |result, comparison| {
            result.and_then(|_| writeln!(f, "{}", comparison))
        })
    }
}

impl Deref for Comparisons {
    type Target = Vec<Comparison>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Comparisons {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::comparisons::tests::common;
    use super::super::UNDEFINED_LABEL_PREFIX;

    fn factorial(n: u32 ) -> u32 {
        match n {
            0   => 1,
            1.. => (1..n+1).product(),
        }
    }
    fn permutations_sample(n: u32) -> u32 {
        match n {
            0|1|2 => factorial(n) - 1,
            3..   => (factorial(n)/(factorial(n-2)*2)),
        }
    }

    fn mock_comparisons(ind_set: &Vec<usize>, allow_self_comparison: bool) -> Comparisons {
        let min_depths = vec![2];
        let names = vec![];
        let genome = Genome::default();
        Comparisons::parse(&ind_set, &min_depths, &names, allow_self_comparison, &genome, 50_000_000)
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

        empty_comp.push(common::mock_comparison(false));
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
    fn display() {
        use std::fmt::Write;
        use super::super::{PAIRS_FORMAT_LEN, COUNT_FORMAT_LEN, AVERG_FORMAT_LEN, DISPL_SEP, FLOAT_FORMAT_PRECISION};

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
                {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}{DISPL_SEP}\
                {: <AVERG_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$}",
                expected_pair_name, 0, 0.0, f64::NAN, f64::NAN, 0.00000, f64::NAN
            ).unwrap();
        }

        assert_eq!(format!("{comparisons}"), expected_output);
    }

}