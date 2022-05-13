use genome::JackknifeBlocks;
use genome::{SNPCoord, Genome};
use crate::pileup::{Pileup, Line, Nucleotide};
use std::{
    fmt,
    collections::BTreeSet
};
use itertools::Itertools;

/// Represents a requested individual Within the pileup.
///  - name      : Name of the individual. Either given through user-input, or constructed as `Ind{index}` by default
///                when no name has been provided.
///  - index     : 0 based index of the individual within the pileup. Note that this index is technically offset by 3,
///                since the first three columns of a pileup respectively define 'chr', 'pos', 'ref'.
///  - min_depth : minimum sequencing depth that is allowed before making a comparison. User-defined, or defaults to 1.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Individual {
    pub name     : String,
    pub index    : usize,
    pub min_depth: u16,
}

impl Individual {
    pub fn new(name: Option<&String>, index: usize, min_depth: u16) -> Individual {
        let name = match name {                            // Parse name and give default if non-existent
            Some(name) => name.to_string(),
            None => format!("Ind{}", index),
        };
        Individual{name, index, min_depth }
    }
}



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
    pair            : (Individual, Individual),
    self_comparison : bool,
    pwd             : u32,
    sum_phred       : u32,
    pub blocks      : JackknifeBlocks,
    pub positions   : BTreeSet<SNPCoord>
}

impl Comparison {
    pub fn new(pair: (Individual, Individual), self_comparison: bool, genome: &Genome, blocksize: u32) -> Comparison {
        Comparison {pair, self_comparison, pwd:0, sum_phred:0, blocks: JackknifeBlocks::new(genome, blocksize), positions: BTreeSet::new()}
    }

    pub fn get_pair_indices(&self) -> [usize; 2] {
        [self.pair.0.index, self.pair.1.index]
    }

    // Check if the sequencing depth of a given overlap is over the minimum required sequencing depth for each individual.
    pub fn satisfiable_depth(&self, pileups: &[Pileup]) -> bool {
        pileups[self.pair.0.index].depth >= self.pair.0.min_depth && pileups[self.pair.1.index].depth >= self.pair.1.min_depth
    }

    // Compare our two individuals at the given SNPposition ; increment the appropriate counters after the comparison 
    // has been made.
    pub fn compare(&mut self, line: &Line) {
        self.positions.insert(line.coordinate.clone());
        let random_nucl: Vec<&Nucleotide> = if self.self_comparison {  // Self comparison => Combination without replacement. 
            line.random_sample_self(&self.pair.0.index)                 //   - possible combinations: n!/(n-2)!
        } else {                                                       // Std comparison  => Permutation.
            line.random_sample_pair(&self.pair)                        //  - possible combinations: n_1 * n_2
        };
        self.add_phred(&random_nucl);

        let current_block = self.blocks.find_block(&line.coordinate);
        current_block.add_count();
        if Self::check_pwd(&random_nucl) {
            self.pwd +=1;
            current_block.add_pwd();
        }
    }

    // Increment our `sum_phred` counter. The incremented value is computed as the average of the two sampled nucleotides.
    fn add_phred(&mut self, nuc: &[&Nucleotide]) {
        self.sum_phred += ( (nuc[0].phred+nuc[1].phred)/2 ) as u32;
    }

    // Check if there is a pairwise difference.
    fn check_pwd(nuc: &[&Nucleotide]) -> bool {
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
        format!("{}-{}", &self.pair.0.name, &self.pair.1.name)
    }
}

impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
            "{: <20} - {: <7} - {: <7} - {: <8.5} - {: <8.5}",
            self.get_pair(),
            self.positions.len(),
            self.pwd,
            self.get_avg_pwd(),
            self.get_avg_phred()
        )
    }
}

/// Vector of all the user-requested comparisons
#[derive(Debug)]
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
        comparisons.push(Comparison::new((pair[0].clone(), pair[1].clone()), self_comparison, genome, blocksize));            
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
    ///Return slice of our Vector of comparisons.
    pub fn get(&self) -> &[Comparison] {
        &self.0
    }

    ///Return mutable slice of our Vector of comparisons.
    pub fn get_mut(&mut self) -> &mut [Comparison] {
        &mut self.0
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

#[cfg(test)]
mod tests {
    use super::*;

    fn factorial(n: i32 ) -> i32 {
        match n {
            0 => 1,
            1.. => (1..n+1).product(),
            _ => panic!("Not an unsigned number")
        }
    }
    fn permutations_sample_2(n: i32) -> i32 {
        match n {
            0.. => (factorial(n)/(factorial((n-2).abs())*2)),
            _ => panic!()
        }
    }

    #[test]
    fn comparisons_length_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let min_depths = vec![2];
            let names = vec![];
            let genome = Genome::default();
            let comparisons = Comparisons::parse(&ind_set, &min_depths, &names, true, &genome, 50_000_000);
            let len = ind_set.len() as i32;
            println!("{:?}", comparisons);
            assert_eq!(comparisons.len() as i32, permutations_sample_2(len)+len); 
        }
    }

    #[test]
    fn comparisons_length_no_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let min_depths = vec![2];
            let names = vec![];
            let genome = Genome::default();
            let comparisons = Comparisons::parse(&ind_set, &min_depths, &names, false, &genome, 50_000_000);
            let len = ind_set.len() as i32;
            assert_eq!(comparisons.len() as i32, permutations_sample_2(len)); 
        }
    }
}