use genome::{Nucleotide, SNPCoord, coordinate::{ChrIdx, Position, Coordinate, derive::Coord}, snp::Allele};
use crate::{comparisons::Individual};
use rand::seq::SliceRandom;

use located_error::prelude::*;

use super::Pileup;

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
#[derive(Debug, Coord)]
pub struct Line {
    pub coordinate : Coordinate,
    pub reference  : Allele,
    pub individuals: Vec<Pileup>,
}

impl Line {
    /// Instantiate a new pileup `Line`
    /// 
    /// # Errors
    /// - [`ParseChr`]   if chromosome is an invalid u8               (field [0] of the pileup)
    /// - [`ParsePos`]   if position   is an invalid u32              (field [1])
    /// - [`ParseRef`]   if reference  is an invalid allele character (field [2])
   /// -  [`ParseDepth`] if any depth  field is an invalid u16        (fields[i%3])
    /// - [`PileupError::RefSkip`] if the pileup line contains reference skips ('[<>]' characters)
    /// - [`PileupError::UnequalLength`] if the base and scores strings do not match in length.
    /// - If any indel is encountered and the program fails to skip it.
    pub fn new(line: &str, ignore_dels: bool) -> Result<Line> {
        use super::PileupError::{ParseChr, ParsePos, ParseRef, ParseDepth};
        let err_context = "While parsing new pileup line";
        let fields: Vec<&str>    = line.split('\t').collect();
        let chromosome: ChrIdx   = fields[0].parse().map_err(ParseChr).loc(err_context)?;
        let position  : Position = fields[1].parse().map_err(ParsePos).loc(err_context)?;
        let reference : Allele   = fields[2].parse().map_err(ParseRef).loc(err_context)?;
        //Loop along individuals
        let mut individuals: Vec<Pileup> = Vec::new();
        for i in (3..fields.len()).step_by(3) {
            let depth           = fields[i].parse().map_err(ParseDepth).loc(err_context)?;
            let (bases, scores) = (fields[i+1], fields[i+2]);

            individuals.push(Pileup::new(reference, depth, bases, scores, ignore_dels)?);
        }
        Ok(Line { coordinate: Coordinate::new(chromosome, position), reference, individuals })
    }

    /// Apply base quality filtering on each individual pileup, according to a given treshold
    /// See: `Pileup::filter_base_quality()`
    pub fn filter_base_quality(&mut self, phred_treshold: u8) {
        self.individuals.iter_mut().for_each(|ind| ind.filter_base_quality(phred_treshold.into()) )
    }

    /// Apply `known_variant` filtering on each individual pileup, according to a given treshold.
    /// See: `Pileup::filter_known_variant()`
    /// 
    /// # Errors
    /// - will bubble any `MissingAltRef` error raised when filtering known variants on a given individual.
    pub fn filter_known_variants(&mut self, known_variant: &SNPCoord) -> Result<()> {
        for individual in &mut self.individuals {
            individual.filter_known_variants(known_variant)
                .with_loc(|| "While filtering known variants" )?;
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
        // @TODO: change to fastrand.
        let mut rng = &mut rand::thread_rng();
        Some(vec![ // @TODO: Yuk! change this horrid allocation.
            self.individuals[pair[0].index].nucleotides.choose(&mut rng)?,
            self.individuals[pair[1].index].nucleotides.choose(&mut rng)?,
        ])
    }

}

#[cfg(test)]
mod tests {
    use anyhow::Result;

    use crate::{pileup};
    use genome::SNPCoord;

    #[test]
    fn line_filter_base_qualities() -> Result<()>{
        println!("Testing Line base_quality filtering.");
        let raw_line="22\t51057923\tC\t6\tTTTTtt\tJEJEEE\t0\t*\t*\t1\tT\tJ";
        let mut line = pileup::Line::new(raw_line, true)?;

        line.filter_base_quality(30);
        assert_eq!(line.individuals[0].get_nucleotides(), "TTTTTT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JEJEEE");

        line.filter_base_quality(40);
        assert_eq!(line.individuals[0].get_nucleotides(), "TT");
        assert_eq!(line.individuals[0].get_scores_ascii(), "JJ");
        Ok(())
    }

    #[test]
    fn line_filter_known_variant_1() -> Result<()> {
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
    fn line_filter_known_variant_2() -> Result<()> {
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