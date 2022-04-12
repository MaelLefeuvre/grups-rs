use crate::genome::{SNPCoord, Chromosome};
use crate::jackknife::JackknifeBlocks;
use std::iter::Peekable;
use std::error::Error;
use rand::seq::SliceRandom;
use std::fmt;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Individual {
    pub name     : String,
    pub index    : usize,
    pub min_depth: u16,
}

impl Individual {
    pub fn new(name: Option<&String>, index: &usize, min_depth: &u16) -> Individual {
        let name = match name {                            // Parse name and give default if non-existent
            Some(name) => name.to_string(),
            None => format!("Ind{}", &index.to_string()),
        };
        Individual {name, index: *index, min_depth: *min_depth }
    }
}

#[derive(Debug)]
pub struct Comparison {
    pub pair           : (Individual, Individual),
    pub self_comparison: bool,
    pub overlap        : u32,
    pub pwd            : u32,
    pub sum_phred      : u32,
    pub blocks         : JackknifeBlocks
}

impl Comparison {
    pub fn new(pair: (Individual, Individual), self_comparison: bool, genome: &Vec<Chromosome>, blocksize: u32) -> Comparison {
        Comparison {pair, self_comparison, overlap: 0, pwd:0, sum_phred:0, blocks: JackknifeBlocks::new(&genome, blocksize)}
    }
    pub fn satisfiable_depth(&self, pileups: &Vec<Pileup>) -> bool {
        pileups[self.pair.0.index].depth >= self.pair.0.min_depth && pileups[self.pair.1.index].depth >= self.pair.1.min_depth
    }

    pub fn compare(&mut self, line: &Line) -> () {
        self.overlap +=1;
        let random_nucl: Vec<&Nucleotide> = if self.self_comparison {                  // Self comparison => Combination without replacement. 
            line.random_sample_self(&self.pair.0.index)                                //   - possible combinations: n!/(n-2)!
        } else {                                                                       // Std comparison  => Permutation.
            line.random_sample_pair(&self.pair)                                        //  - possible combinations: n_1 * n_2
        };
        self.add_phred(&random_nucl);

        let current_block = self.blocks.find_block(&line.coordinate);
        current_block.add_count();
        if Self::check_pwd(&random_nucl) {
            self.pwd +=1;
            current_block.add_pwd();
        }
        ()
    }

    fn add_phred(&mut self, nuc: &Vec<&Nucleotide>) -> () {
        self.sum_phred += ( (nuc[0].phred+nuc[1].phred)/2 ) as u32;
    }

    fn check_pwd(nuc: &Vec<&Nucleotide>) -> bool {
        nuc[0].base != nuc[1].base
    }

    pub fn get_avg_pwd(&self) -> f64 {
        self.pwd as f64 / self.overlap as f64
    }

    pub fn get_avg_phred(&self) -> f64 {
        self.sum_phred as f64/self.overlap as f64
    }

    pub fn get_pair (&self,)-> String {
        format!("{}-{}", self.pair.0.name, self.pair.1.name).to_owned()
    }
}

impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
            "{: <20} - {: <7} - {: <7} - {: <8.5} - {: <8.5}",
            self.get_pair(),
            self.overlap,
            self.pwd,
            self.get_avg_pwd(),
            self.get_avg_phred()
        )
    }
}

/// Simple struct representing a given nucleotides.
/// BQ scores are housed in phred-33 scale format.
#[derive(Debug)]
pub struct Nucleotide {
    pub base: char,
    pub phred: u8,
}

impl Nucleotide {
    fn new(base: &char, score: &char) -> Nucleotide {
        let phred = Nucleotide::to_phred(score);            
        Nucleotide {base: *base, phred: phred}
    }

    /// Convert the BQ score back to the ASCII format
    pub fn get_score_ascii(&self) -> char {
        (self.phred + 33) as char 
    }

    /// Convert an ASCII BQ score to phred-33.
    /// Used during construction.
    fn to_phred(score: &char) -> u8 {
        let phred_scale: u8 =33;
        (*score as u8) - phred_scale
    }
}


#[derive(Debug)]
pub enum PileupError {
    RefSkipError,
    LengthError,
    ParseLineError,
}

impl Error for PileupError {}

impl fmt::Display for PileupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::RefSkipError   => write!(f, "Cannot handle reference skips within pileup. ('>' '<')"),
            Self::LengthError    => write!(f, "Length of nucleotides and scores differ vector differ."),
            Self::ParseLineError => write!(f, "Failed to parse pileup line"),
        }
    }
}


/// Pileup record of a single individual, at a given position.
/// Nested structure: Pileup +-> depth
///                          L-> Vec<Nucleotides> +-> base
///                                               L-> score
#[derive(Debug)]
pub struct Pileup {
    pub depth: u16,
    pub nucleotides: Vec<Nucleotide>,
}

impl Pileup {
    pub fn new(depth: u16, bases: String, scores: String, ignore_dels: bool) -> Result<Pileup, Box<dyn Error>> {

        let bases = &bases.to_uppercase();    // Convert antisense to forward
        let bases = &bases.replace(",", "."); //

        // Loop along nucleotides.
        let mut nucleotides: Vec<Nucleotide> = Vec::new();
        let mut scores_vec = scores.chars();

        // Filter out non-selectible nucleotides.
        let mut chars = bases.chars().peekable();
        while let Some(n) = chars.by_ref().next() {
            match n {
                '+'|'-' => {Self::skip_indel(&mut chars); continue},        //Skip indels
                '^'     => {chars.nth(0); continue},                        //Skip starts
                '$'     => {continue},                                      //Skip end
                '*'     => if ignore_dels {continue},                       //Skip deletion if required
                '>'|'<' => return Err(Box::new(PileupError::RefSkipError)), //Bail if we found a refskip
                _       => ()
            }
            match &scores_vec.next() {
                Some(score) => nucleotides.push(Nucleotide::new(&n, &score)),
                None => return Err(Box::new(PileupError::LengthError))
            };
        }
        Ok(Pileup { depth, nucleotides })
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
    /// TODO : Better error handling for numeric checking
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
///        use .nth(0) instead on chrom, pos, ref, 
///        instead of collecting everything and parsing individually
///        => Except for 'reference' we should return a ParseError i
///           f we encounter None at any point
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
            let depth :    u16 = split_line[i].parse()?;
            let bases : String = split_line[i+1].parse()?;
            let scores: String = split_line[i+2].parse()?;

            individuals.push(Pileup::new(depth, bases, scores, ignore_dels)?);
        }
        Ok(Line {
            coordinate: SNPCoord{chromosome, position, reference, alternate: None },
            individuals: individuals
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
    /// Two nucleotide sampled, without replacement.
    pub fn random_sample_self(&self, index: &usize) -> Vec<&Nucleotide> {
        let mut rng = &mut rand::thread_rng();
        self.individuals[*index].nucleotides.choose_multiple(&mut rng, 2).collect() 
    }

    /// Apply random sampling on a a pair of individuals for pairwise-comparison.
    /// One nucleotide sampled per-ind, without replacement.
    pub fn random_sample_pair(&self, pair: &(Individual, Individual)) -> Vec<&Nucleotide> {
        let mut rng = &mut rand::thread_rng();
        vec![
            self.individuals[pair.0.index].nucleotides.choose(&mut rng).unwrap(),
            self.individuals[pair.1.index].nucleotides.choose(&mut rng).unwrap(),
        ]

    }

}

#[cfg(test)]
mod tests {
    use crate::pileup;
    use crate::genome::SNPCoord;

    fn create_dummy_pileup(line_input: &str, score_input: &str, ignore_dels: bool) -> pileup::Pileup {
        let line_input : String = line_input.to_string();
        let scores_input : String = score_input.to_string();
        return pileup::Pileup::new(line_input.len() as u16, line_input, scores_input, ignore_dels).unwrap();
    }

    #[test]
    fn line_filter_base_qualities() {
        println!("Testing Line base_quality filtering.");
        let raw_line="22\t51057923\tC\t6\tTTTTtt\tJEJEEE\t0\t*\t*\t1\tT\tJ";
        let mut line = pileup::Line::new(&raw_line, true).unwrap();
        // {
        //    Ok(line) => line,
        //    Err(e)   => {log::error!("Error: {e}", e)},
        //};

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


    #[test]
    fn pileup_mutate_for_rev_reference() {
        println!("Testing reverse to forward conversion: ',' -> '.'");
        let line_input   = "...,..,.,...,.";
        let line_expect  = "..............";
        let scores_input = "JEJEEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false);
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_mutate_for_rev_alternate() {
        println!("Testing reverse to forward conversion: [atcgn] -> [ATCGN] ");
        let line_input   = ",.actgn,.NGTCA,.";
        let line_expect  = "..ACTGN..NGTCA..";
        let scores_input = "JEJEEECc$cacJgGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false);
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_start_stop() {
        println!("Testing start/stop filtering");
        let line_input   = ",.$ac.N^JTA$^AC,.";
        let line_expect  = "..AC.NTAC..";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, false);
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_missing() {
        println!("Testing missing base filtering");
        let line_input   = ",.$a*c.N^JTA*C,.";
        let line_expect  = "..AC.NTAC..";
        let scores_input = "JEECc$cagGg";
        let pileup = create_dummy_pileup(line_input, scores_input, true);
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }

    #[test]
    fn pileup_filter_indels() {
        println!("Testing indel filtering");
        let line_input   = ",..,,+4ACTAGca,,.,-2AT..,.+15ATCGCCCCGCCCTAGc";
        let line_expect  = ".....GCA........C";
        let scores_input = "JEEeCCeCCc$cagGgc";
        let pileup = create_dummy_pileup(line_input, scores_input, true);
        println!("{:?}\n{:?}", pileup.get_nucleotides(), line_expect);
        assert_eq!(pileup.get_nucleotides(), line_expect);
    }
}
