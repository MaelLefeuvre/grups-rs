use genome::{SNPCoord, coordinate::{ChrIdx, Position, Coordinate}, snp::Allele};
use crate::{comparisons::Individual, io::SNPReaderError};
use std::error::Error;
use rand::seq::SliceRandom;

use located_error::prelude::*;

use super::{Nucleotide, Pileup};

/// Parsed line of a pileup file entry, containing the both coordinates and 
/// Pileup of each individual.
/// 
/// This structure is heavily nested:
/// `Line` +-> `SNPCoord` +-> chromosome
///        |              +-> position
///        |              +-> ref
///        |              L-> alt
///        L-> `Vec<Pileup>` +-> depth
///                          L-> `Vec<Nucleotides>` +-> base
///                                                 L-> score
/// 
#[derive(Debug)]
pub struct Line {
    pub coordinate : Coordinate,
    pub reference  : Allele,
    pub individuals: Vec<Pileup>,
}

impl Line {
    /// Instantiate a new pileup `Line`
    /// 
    /// # Errors
    /// - `ParseIntError` if chromosome and position fails to get parsed (fields [0] and [1] of the pileup)
    /// - `ParseIntError` if any of the `depth` fields fails to get parsed into an integer
    ///    these fields are located at the indices where `(i+3) %% 3 == 0`
    /// - `PileupError::RefSkip` if the pileup line contains reference skips ('[<>]' characters)
    /// - `PileupError::UnequalLength` if the base and scores strings do not match in length.
    /// - If any indel is encountered and the program fails to skip it.
    pub fn new(line: &str, ignore_dels: bool) -> Result<Line> {
        let split_line: Vec<&str>    = line.split('\t').collect();
        let chromosome: ChrIdx       = split_line[0].parse()?;
        let position  : Position     = split_line[1].parse()?;
        let reference : Allele       = split_line[2].parse().loc(format!("While parsing reference allele {}", split_line[2]))?;
        //Loop along individuals
        let mut individuals: Vec<Pileup> = Vec::new();
        for i in (3..split_line.len()).step_by(3) {
            let depth : u16  = split_line[i].parse()?;
            let bases : &str = split_line[i+1];
            let scores: &str = split_line[i+2];

            individuals.push(Pileup::new(reference, depth, bases, scores, ignore_dels)?);
        }
        Ok(Line {
            coordinate: Coordinate::new(chromosome, position),
            reference,
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

    /// Apply `known_variant` filtering on each individual pileup, according to a given treshold.
    /// See: `Pileup::filter_known_variant()`
    /// 
    /// # Errors
    /// - will bubble any `MissingAltRef` error raised when filtering known variants on a given individual.
    pub fn filter_known_variants(&mut self, known_variant: &SNPCoord) -> Result<()> {
        for individual in &mut self.individuals {
            individual.filter_known_variants(known_variant)?;
        }
        Ok(())
    }

    /// Apply random sampling on a single individual for self-comparison.
    /// Two nucleotide sampled (without replacement).
    #[must_use]
    pub fn random_sample_self(&self, index: &usize) -> Vec<&Nucleotide> {
        let mut rng = &mut rand::thread_rng();
        self.individuals[*index].nucleotides.choose_multiple(&mut rng, 2).collect() 
    }

    /// Apply random sampling on a a pair of individuals for pairwise-comparison.
    /// One nucleotide sampled per individual (without replacement).
    /// @TODO: rng should not be initialized within function.
    #[must_use]
    pub fn random_sample_pair(&self, pair: &[Individual; 2]) -> Option<Vec<&Nucleotide>> {
        let mut rng = &mut rand::thread_rng();
        Some(vec![
            self.individuals[pair[0].index].nucleotides.choose(&mut rng)?,
            self.individuals[pair[1].index].nucleotides.choose(&mut rng)?,
        ])
    }

}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use crate::{pileup};
    use genome::SNPCoord;

    #[test]
    fn line_filter_base_qualities() -> Result<(), Box<dyn Error>>{
        println!("Testing Line base_quality filtering.");
        let raw_line="22\t51057923\tC\t6\tTTTTtt\tJEJEEE\t0\t*\t*\t1\tT\tJ";
        let mut line = pileup::Line::new(raw_line, true)?;

        line.filter_base_quality(&30);
        assert_eq!(line.individuals[0].get_nucleotides(), "TTTTTT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JEJEEE");

        line.filter_base_quality(&40);
        assert_eq!(line.individuals[0].get_nucleotides(), "TT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JJ");
        Ok(())
    }

    #[test]
    fn line_filter_known_variant_1() -> Result<(), Box<dyn Error>> {
        println!("Testing Line known_variant filtration.");
        let raw_line="2\t21303470\tN\t0\t*\t*\t8\tTTcTTtt^Ft\tEEJEEEEE\t0\t*\t*";
        let mut line = pileup::Line::new(raw_line, true)?;

        let known_variant = SNPCoord::try_new(2, 21303470, 'C', 'A')?;
        line.filter_known_variants(&known_variant)?;
        assert_eq!(line.individuals[1].get_nucleotides(), String::from('C'));
        assert_eq!(line.individuals[1].get_scores_ascii(), String::from('J'));
        assert_eq!(line.individuals[1].depth, 1);
        Ok(())
    }

    #[test]
    fn line_filter_known_variant_2() -> Result<(), Box<dyn Error>> {
        println!("Testing Line known_variant filtration.");
        let raw_line="2\t21303470\tT\t0\t*\t*\t8\t..c..,,^F,\tEEJEEEEE\t0\t*\t*";
        let mut line = pileup::Line::new(raw_line, true)?;

        let known_variant = SNPCoord::try_new(2, 21303470, 'T', 'A')?;
        line.filter_known_variants(&known_variant)?;
        assert_eq!(line.individuals[1].get_nucleotides(),  String::from("TTTTTTT"));
        assert_eq!(line.individuals[1].get_scores_ascii(), String::from("EEEEEEE"));
        assert_eq!(line.individuals[1].depth, 7);
        Ok(())
    }
}