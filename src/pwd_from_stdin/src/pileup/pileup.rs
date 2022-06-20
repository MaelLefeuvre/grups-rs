use genome::SNPCoord;
use std::ops::Add;
use std::{iter::Peekable, collections::HashMap};
use std::error::Error;
use itertools::Itertools;

use super::{Nucleotide, PileupError};




/// Pileup record of a single individual, at a given position.
/// Nested structure: Pileup +-> depth
///                          L-> Vec<Nucleotides> +-> base
///                                               L-> score
/// 
/// # TODO: migrate nucleotide formating from `Pileup::new()` to `Nucleotide::new()`
/// this will greatly increase readability and dilute responsibility.
#[derive(Debug)]
pub struct Pileup {
    pub depth: u16,
    pub nucleotides: Vec<Nucleotide>,
}

impl Pileup {
    pub fn new(depth: u16, bases: &str, scores: &str, ignore_dels: bool) -> Result<Pileup, Box<dyn Error>> {

        let bases = &bases.to_uppercase();    // Convert antisense to forward
        let bases = &bases.replace(',', "."); //

        // Loop along nucleotides.
        let mut nucleotides: Vec<Nucleotide> = Vec::new();
        let mut scores_vec = scores.chars();

        // Filter out non-selectible nucleotides.
        let mut chars = bases.chars().peekable();
        while let Some(n) = chars.by_ref().next() {
            match n {
                '+'|'-' => {Self::skip_indel(&mut chars); continue},   //Skip indels
                '^'     => {chars.next(); continue},                   //Skip starts
                '$'     => {continue},                                 //Skip end
                '*'     => if ignore_dels {continue},                  //Skip deletion if required
                '>'|'<' => return Err(Box::new(PileupError::RefSkip)), //Bail if we found a refskip
                _       => ()
            }
            match scores_vec.next() {
                Some(score) => nucleotides.push(Nucleotide::new(n, score)),
                None => return Err(Box::new(PileupError::UnequalLength))
            };
        }
        Ok(Pileup { depth, nucleotides })
    }

    pub fn observation_set(&self) -> (HashMap<char, f64>, f64) {
        let mut set: HashMap<char, f64> = HashMap::new();
        let mut sum_phred = 0;
        for nucleotide in self.nucleotides.iter() {
            *set.entry(nucleotide.base).or_insert(0.0) += 1.0;
            sum_phred += nucleotide.phred;
        }
        for count in set.values_mut() {
            *count = *count / self.nucleotides.len() as f64;
        }
        let avg_phred = sum_phred as f64 / self.nucleotides.len() as f64 ;
        (set, avg_phred)
    }

    /// Apply base quality filtration for each Nucleotide, using a given a phred treshold
    pub fn filter_base_quality(&mut self, phred_treshold: &u8) {
        self.nucleotides.retain(|nucleotide| {
            nucleotide.phred >= *phred_treshold
        });
        self.update_depth();
    }

    /// Apply triallelic filtration for each Nucleotide, using a given SNPCoordinate.
    /// Note that SNPCoordinates must contain values for their `reference` and `alternate` fields
    /// 
    /// TODO : Better error handling for known_variant unwrapping.
    ///        => Do .ok() -> early return if None
    pub fn filter_known_variants(&mut self, known_variant: &SNPCoord) {
        let (reference, alternate) = (known_variant.reference.unwrap(), known_variant.alternate.unwrap());
        self.nucleotides.retain(|nucleotide| {
            vec!['.', reference, alternate].contains(&nucleotide.base)
        });
        self.update_depth();
    }

    /// Return a string of nucleotides. Mainly for Tests and Debug.
    pub fn get_nucleotides(&self) -> String {
        let mut pileup_nuc = String::from("");
        for nuc in &self.nucleotides {
            pileup_nuc.push(nuc.base);
        }
        pileup_nuc
    }

    /// Return a string of BQ scores in ASCII format. Mainly for Tests and Debug.
    #[must_use]
    pub fn get_scores_ascii(&self) -> String {
        let mut pileup_nuc = String::from("");
        for nuc in &self.nucleotides {
            pileup_nuc.push(nuc.get_score_ascii());
        }
        pileup_nuc
    }

    /// Run through a Peekable Iterator of nucleotides to parse the number of characters that should be skipped
    /// When encountering an indel.
    /// Indel regex: [+-][0-9]+[ATCGANatcgan]+
    ///               --  ---   ------------
    ///               |   |     + Sequence
    ///               |   + Length of Sequence
    ///               + Identifier ('-' = deletion , '+' = insertion=
    /// 
    /// TODO : Better error handling for numeric checking. Remove those nasty `.unwrap()`
    fn skip_indel<I: Iterator<Item = char>>(chars: &mut Peekable<I>){
        let mut digit: String = "".to_string();
        while chars.peek().unwrap().is_numeric() {
            digit.push(chars.next_if(|&x| x.is_numeric()).unwrap());
        }
        let skip = digit.parse::<usize>().unwrap() -1 ;
        chars.nth(skip);
    }

    /// Refresh the depth after quality filtration.
    /// Mainly used by: - self.filter_known_variants()
    ///                 - self.filter_base_quality()
    fn update_depth(&mut self) {
        self.depth = self.nucleotides.len() as u16;

    }
}

#[cfg(test)]
mod tests {

    use super::*;

    fn create_dummy_pileup(line_input: &str, score_input: &str, ignore_dels: bool) -> Result<Pileup, Box<dyn Error>> {
        return Pileup::new(line_input.len() as u16, line_input, score_input, ignore_dels);
    }

    #[test]
    fn pileup_mutate_for_rev_reference() {
        println!("Testing reverse to forward conversion: ',' -> '.'");
        let line_input   = "...,..,.,...,.";
        let line_expect  = "..............";
        let scores_input = "JEJEEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false).unwrap();
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_mutate_for_rev_alternate() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.actgn,.NGTCA,.";
        let line_expect  = "..ACTGN..NGTCA..";
        let scores_input = "JEJEEECc$cacJgGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false).unwrap();
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_start_stop() {
        println!("Testing start/stop filtering");
        let line_input   = ",.$ac.N^JTA$^AC,.";
        let line_expect  = "..AC.NTAC..";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false).unwrap();
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_missing() {
        println!("Testing missing base filtering");
        let line_input   = ",.$a*c.N^JTA*C,.";
        let line_expect  = "..AC.NTAC..";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, true).unwrap();
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_indels() {
        println!("Testing indel filtering");
        let line_input   = ",..,,+4ACTAGca,,.,-2AT..,.+15ATCGCCCCGCCCTAGc";
        let line_expect  = ".....GCA........C";
        let scores_input = "JEEeCCeCCc$cagGgc";
        let pileup = create_dummy_pileup(line_input, scores_input, true).unwrap();
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_refskip_error() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.act>gn,.NGTCA,.";
        let scores_input = "JEJEEECc$cacJgGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false);
        assert!(pileup.is_err())
    }

    #[test]
    fn pileup_unequal_length_error() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.actgn,.NGTCA,.";
        let scores_input = "JEcacJGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false);
        assert!(pileup.is_err())
    }

}