use genome::SNPCoord;
use crate::comparisons::Individual;
use std::error::Error;
use rand::seq::SliceRandom;

use super::{Nucleotide, Pileup};

/// Parsed line of a pileup file entry, containing the both coordinates and 
/// Pileup of each individual.
/// 
/// This structure is heavily nested:
/// Line +-> SNPCoord +-> chromosome
///      |            +-> position
///      |            +-> ref
///      |            L-> alt
///      L-> Vec<Pileups> +-> depth
///                       L-> Vec<Nucleotides> +-> base
///                                            L-> score
/// 
/// TODO : Error handling should be better taken of. 
///        use .nth(0) or .next() instead on chrom, pos, ref, 
///        instead of collecting everything and parsing individually
///        => Except for 'reference' we should return a ParseError i
///           f we encounter None at any point
/// 
#[derive(Debug)]
pub struct Line {
    pub coordinate : SNPCoord,
    pub individuals: Vec<Pileup>,
}

impl Line {
    pub fn new(line: &str, ignore_dels: bool) -> Result<Line, Box<dyn Error>> {
        let split_line: Vec<&str>    = line.split('\t').collect();
        let chromosome: u8           = split_line[0].parse()?;
        let position  : u32          = split_line[1].parse()?;
        let reference : Option<char> = split_line[2].parse().ok();
        //Loop along individuals
        let mut individuals: Vec<Pileup> = Vec::new();
        for i in (3..split_line.len()).step_by(3) {
            let depth : u16  = split_line[i].parse()?;
            let bases : &str = split_line[i+1];
            let scores: &str = split_line[i+2];

            individuals.push(Pileup::new(depth, bases, scores, ignore_dels)?);
        }
        Ok(Line {
            coordinate: SNPCoord{chromosome, position, reference, alternate: None },
            individuals
        })
    }

    /// Apply base quality filtering on each individual pileup, according to a given treshold
    /// See: `Pileup::filter_base_quality()`
    pub fn filter_base_quality(&mut self, phred_treshold: &u8) {
        for individual in &mut self.individuals {
            individual.filter_base_quality(phred_treshold);
        }
    }

    /// Apply known_variant filtering on each individual pileup, according to a given treshold.
    /// See: Pileup::filter_known_variant()
    pub fn filter_known_variants(&mut self, known_variant: &SNPCoord) {
        for individual in &mut self.individuals {
            individual.filter_known_variants(known_variant);
        }
    }

    /// Apply random sampling on a single individual for self-comparison.
    /// Two nucleotide sampled (without replacement).
    pub fn random_sample_self(&self, index: &usize) -> Vec<&Nucleotide> {
        let mut rng = &mut rand::thread_rng();
        self.individuals[*index].nucleotides.choose_multiple(&mut rng, 2).collect() 
    }

    /// Apply random sampling on a a pair of individuals for pairwise-comparison.
    /// One nucleotide sampled per individual (without replacement).
    pub fn random_sample_pair(&self, pair: &[Individual; 2]) -> Vec<&Nucleotide> {
        let mut rng = &mut rand::thread_rng();
        vec![
            self.individuals[pair[0].index].nucleotides.choose(&mut rng).unwrap(),
            self.individuals[pair[1].index].nucleotides.choose(&mut rng).unwrap(),
        ]

    }

}

#[cfg(test)]
mod tests {
    use crate::pileup;
    use genome::SNPCoord;

    #[test]
    fn line_filter_base_qualities() {
        println!("Testing Line base_quality filtering.");
        let raw_line="22\t51057923\tC\t6\tTTTTtt\tJEJEEE\t0\t*\t*\t1\tT\tJ";
        let mut line = pileup::Line::new(&raw_line, true).unwrap();

        line.filter_base_quality(&30);
        assert_eq!(line.individuals[0].get_nucleotides(), "TTTTTT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JEJEEE");

        line.filter_base_quality(&40);
        assert_eq!(line.individuals[0].get_nucleotides(), "TT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JJ");
    }

    #[test]
    fn line_filter_known_variant_1() {
        println!("Testing Line known_variant filtration.");
        let raw_line="2\t21303470\tN\t0\t*\t*\t8\tTTcTTtt^Ft\tEEJEEEEE\t0\t*\t*";
        let mut line = pileup::Line::new(&raw_line, true).unwrap();

        let known_variant = SNPCoord { chromosome: 2, position: 21303470, reference: Some('C'), alternate: Some('A') };
        line.filter_known_variants(&known_variant);
        assert_eq!(line.individuals[1].get_nucleotides(), String::from('C'));
        assert_eq!(line.individuals[1].get_scores_ascii(), String::from('J'));
        assert_eq!(line.individuals[1].depth, 1);
    }

    #[test]
    fn line_filter_known_variant_2() {
        println!("Testing Line known_variant filtration.");
        let raw_line="2\t21303470\tT\t0\t*\t*\t8\t..c..,,^F,\tEEJEEEEE\t0\t*\t*";
        let mut line = pileup::Line::new(&raw_line, true).unwrap();

        let known_variant = SNPCoord { chromosome: 2, position: 21303470, reference: Some('T'), alternate: Some('A') };
        line.filter_known_variants(&known_variant);
        assert_eq!(line.individuals[1].get_nucleotides(),  String::from("......."));
        assert_eq!(line.individuals[1].get_scores_ascii(), String::from("EEEEEEE"));
        assert_eq!(line.individuals[1].depth, 7);
    }
}