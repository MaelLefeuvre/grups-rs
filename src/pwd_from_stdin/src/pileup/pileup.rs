use genome::{SNPCoord, snp::Allele};
use std::{iter::Peekable, collections::HashMap, str::FromStr};
//use std::error::Error;
use located_error::prelude::*;
use crate::io::SNPReaderError;

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
    pub nucleotides: Vec<Nucleotide>
}

impl Pileup {
    /// Instantiate a new `Pileup` entry, tied to a given `Individual`.
    /// 
    /// # Errors 
    /// - If an indel is encountered and the program fails to skip it.
    /// -  `RefSkip` if a reference skip is encountered ('[<>]' characters).
    /// - `UnequalLength` if the base and scores strings do not match in length.
    pub fn new(reference: Allele, depth: u16, bases: &str, scores: &str, ignore_dels: bool) -> Result<Pileup> {

        let bases = &bases.to_uppercase();    // Convert antisense to forward
        //let bases = &bases.replace(',', "."); //

        // Loop along nucleotides.
        let mut nucleotides = Vec::with_capacity(scores.len());
        let mut scores_vec  = scores.chars();

        // Filter out non-selectible nucleotides.
        let mut chars = bases.chars().peekable();
        while let Some(n) = chars.by_ref().next() {
            match n {
                '+'|'-' => {Self::skip_indel(&mut chars)?; continue},  // Skip indels
                '^'     => {chars.next(); continue},                   // Skip starts
                '$'     => {continue},                                 // Skip end
                '*'     => if ignore_dels {continue},                  // Skip deletion if required
                '>'|'<' => return Err(anyhow!(PileupError::RefSkip)),  //Bail if we found a refskip
                _       => ()
            }
            match scores_vec.next() {
                Some(score) => {
                    let n = match n {
                        ',' | '.' => Nucleotide::new(reference, score), // These symbols imply it's the reference allele
                        _         => Nucleotide::try_new(n, score).with_loc(||format!("While parsing {n}"))?,          // Else, try parsing the allele from the given char.
                    };
                    nucleotides.push(n)
                },
                None => return Err(anyhow!(PileupError::UnequalLength))
            };
        }
        Ok(Pileup { depth, nucleotides })
    }

    /// Generate a table of frequency of nucleotide observations for this pileup entry.
    /// Returns:
    ///  - a frequency table as a `HashMap` (K: nucleotide, V: frequency of obs)
    ///  - the `avg_phred` score for all observed nucleotides.
    #[must_use]
    pub fn observation_set(&self) -> (HashMap<Allele, f64>, f64) {
        let mut set: HashMap<Allele, f64> = HashMap::new();
        let mut sum_phred: u32 = 0;
        for nucleotide in &self.nucleotides {
            *set.entry(nucleotide.base).or_insert(0.0) += 1.0;
            sum_phred += u32::from(nucleotide.phred);
        }
        for count in set.values_mut() {
            *count /= self.nucleotides.len() as f64;
        }
        let avg_phred = f64::from(sum_phred) / self.nucleotides.len() as f64 ;
        (set, avg_phred)
    }

    /// Apply base quality filtration for each Nucleotide, using a given a phred treshold
    pub fn filter_base_quality(&mut self, phred_treshold: &u8) {
        self.nucleotides.retain(|nucleotide| {
            nucleotide.phred >= *phred_treshold
        });
        self.update_depth();
    }

    /// Apply triallelic filtration for each Nucleotide, using a given `SNPCoord`.
    /// 
    /// Note that the provided `SNPCoord` must contain values for their `reference` 
    /// and `alternate` fields.
    /// 
    /// # Errors
    ///  - Throws `MissingAltRef` if the `known_variant` does not contain REF/ALT information.
    /// 
    pub fn filter_known_variants(&mut self, known_variant: &SNPCoord) -> Result<()> {
        use crate::io::SNPReaderError::MissingAltRef;
        let (reference, alternate) = (
            known_variant.reference,//.ok_or(MissingAltRef(*known_variant))?,
            known_variant.alternate//.ok_or(MissingAltRef(*known_variant))?
        );
        self.nucleotides.retain(|nucleotide| {
            [reference, alternate].contains(&nucleotide.base)
        });
        self.update_depth();
        Ok(())
    }

    /// Return a string of nucleotides. Mainly for Tests and Debug.
    #[must_use]
    #[cfg(test)]
    pub fn get_nucleotides(&self) -> String {
        let mut pileup_nuc = String::from("");
        for nuc in &self.nucleotides {
            pileup_nuc.push(char::from(&nuc.base));
        }
        pileup_nuc
    }

    /// Return a string of BQ scores in ASCII format. Mainly for Tests and Debug.
    #[must_use]
    #[cfg(test)]
    pub fn get_scores_ascii(&self) -> String {
        let mut pileup_nuc = String::from("");
        for nuc in &self.nucleotides {
            pileup_nuc.push(nuc.get_score_ascii());
        }
        pileup_nuc
    }

    /// Run through a Peekable Iterator of nucleotides to parse the number of characters that should be skipped
    /// When encountering an indel.
    /// 
    /// 
    /// ```Text
    /// Indel regex: `[+-][0-9]+[ATCGANatcgan]+`
    ///                --  ---   ------------
    ///                |   |     | 
    ///                |   |     +- Sequence
    ///                |   | 
    ///                |   +- Length of Sequence
    ///                | 
    ///                +- Identifier ('+' = insertion, '-' = deletion)
    /// ```
    /// # Errors
    /// - If the length of the sequence (following [+-]) cannot get parsed into an integer.
    fn skip_indel<I: Iterator<Item = char>>(chars: &mut Peekable<I>) -> Result<()> {
        let mut digit: String = "".to_string();
        let err_msg = "Error while skipping indels: characters following indel identifier are not numeric.";
        while chars.peek().loc(err_msg)?.is_numeric() {
            digit.push(
                chars.next_if(|&x| x.is_numeric())
                    .loc(err_msg)?
            );
        }

        // Convert digit from string to numeric. Throw an error if this fails.
        let digit = digit.parse::<usize>()
            .with_loc(|| {
                format!("Invalid indel identifier within the pileup file. \
                    Got [{{err}}] while attempting to skip indel"
            )
        })?;
        let skip = digit -1 ;
        chars.nth(skip);
        Ok(())
    }

    /// Refresh the depth after quality filtration.
    /// Mainly used by: - `self.filter_known_variants()`
    ///                 - `self.filter_base_quality()`
    fn update_depth(&mut self) {
        self.depth = self.nucleotides.len() as u16;

    }
}

#[cfg(test)]
mod tests {

    use super::*;

    fn create_dummy_pileup(reference: Allele, line_input: &str, score_input: &str, ignore_dels: bool) -> Result<Pileup> {
        return Pileup::new(reference, line_input.len() as u16, line_input, score_input, ignore_dels);
    }

    #[test]
    fn pileup_mutate_for_rev_reference() -> Result<()> {
        println!("Testing reverse to forward conversion: ',' -> '.'");
        let line_input   = "...,..,.,...,.";
        let line_expect  = "CCCCCCCCCCCCCC";
        let scores_input = "JEJEEECc$cagGg";
        let pileup = create_dummy_pileup(Allele::C, line_input, scores_input, false)?;
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
        Ok(())
    }

    #[test]
    fn pileup_mutate_for_rev_alternate() -> Result<()> {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.actgn,.NGTCA,.";
        let line_expect  = "AAACTGNAANGTCAAA";
        let scores_input = "JEJEEECc$cacJgGg";
        let pileup = create_dummy_pileup(Allele::A, line_input, scores_input, false)?;
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
        Ok(())
    }

    #[test]
    fn pileup_filter_start_stop() -> Result<()> {
        println!("Testing start/stop filtering");
        let line_input   = ",.$ac.N^JTA$^AC,.";
        let line_expect  = "NNACNNTACNN";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(Allele::N, line_input, scores_input, false)?;
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
        Ok(())
    }

    #[test]
    fn pileup_filter_missing() -> Result<()> {
        println!("Testing missing base filtering");
        let line_input   = ",.$a*c.N^JTA*C,.";
        let line_expect  = "NNACNNTACNN";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(Allele::N, line_input, scores_input, true)?;
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
        Ok(())
    }

    #[test]
    fn pileup_filter_indels() -> Result<()> {
        println!("Testing indel filtering");
        let line_input   = ",..,,+4ACTAGca,,.,-2AT..,.+15ATCGCCCCGCCCTAGc";
        let line_expect  = "NNNNNGCANNNNNNNNC";
        let scores_input = "JEEeCCeCCc$cagGgc";
        let pileup = create_dummy_pileup(Allele::N, line_input, scores_input, true)?;
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
        Ok(())
    }

    #[test]
    fn pileup_refskip_error() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.act>gn,.NGTCA,.";
        let scores_input = "JEJEEECc$cacJgGg";
        let pileup = create_dummy_pileup(Allele::N, line_input, scores_input, false);
        assert!(pileup.is_err())
    }

    #[test]
    fn pileup_unequal_length_error() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.actgn,.NGTCA,.";
        let scores_input = "JEcacJGg";
        let pileup = create_dummy_pileup(Allele::N, line_input, scores_input, false);
        assert!(pileup.is_err())
    }

}