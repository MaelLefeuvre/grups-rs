use std::{
    collections::BTreeMap,
    ops::{Deref, DerefMut}, 
    error::Error,
    path::Path,
    io::{BufRead, BufReader},
    fs::File
};

use crate::{
    geneticmap::GeneticMap,
    gamete::Gamete,
    chromosome::Chromosome,
};
use log::{warn, info, debug};

#[derive(Debug)]
pub enum FastaIndexReaderError {
    FileNotFoundError(String, String),
    ParseIntError(String, u8, String),
    InvalidFileFormatError(String)
}
impl Error for FastaIndexReaderError {}

impl std::fmt::Display for FastaIndexReaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::FileNotFoundError(path, err)   => write!(f, "{}: {}",path, err),
            Self::ParseIntError(path, line, err) => write!(f, "{}: Line {} - Invalid line format [{}]", path, line, err),
            Self::InvalidFileFormatError(ext)  => write!(f, "Invalid fasta file format: expected [.fasta|.fa], got '{}'", ext),
        }
    }
}



#[derive(Debug, Clone)]
pub struct Genome(BTreeMap<u8, Chromosome>);

impl Deref for Genome {
    type Target = BTreeMap<u8, Chromosome>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genome {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Genome {

    pub fn new() -> Genome {
        Genome(BTreeMap::new())
    }    

    pub fn from(chromosomes: &[Chromosome]) -> Genome {
        let mut genome = Genome(BTreeMap::new());
        for chr in chromosomes {
            genome.add_chromosome(chr.index, chr.name, chr.length);
        }
        genome
    }

    /// Read a `.fasta.fai` file and parse it into a vector of Chromosome structs.
    /// 
    /// TODO: - convert 'path' variable to &str
    ///       - Check for .fai.fai extension doubles.
    pub fn from_fasta_index(path: &str) -> Result<Genome,FastaIndexReaderError> {
        use FastaIndexReaderError::{InvalidFileFormatError, FileNotFoundError, ParseIntError};
        info!("Parsing reference genome: {}", path);

        let fai = match Path::new(path).extension() { // Check file extension
            None                 => return Err(InvalidFileFormatError("None".to_string())),
            Some(osstr)          => match osstr.to_str(){ 
                Some("fai")          => path.to_string(),                          
                Some("fasta" | "fa") => format!("{}{}", path, ".fai"), // Append '.fai' if it is not there
                Some(other)    => return Err(InvalidFileFormatError(".".to_string()+other)),
                None                 => return Err(InvalidFileFormatError("None".to_string())),
            }        
        };
        debug!("fasta.fai file : {}", fai);
    
    
        let mut genome = Self::new();
        let file = BufReader::new(match File::open(fai.clone()) {
            Ok(file) => file,
            Err(err) => return Err(FileNotFoundError(fai, err.to_string())),
        });
    
        let mut skipped_chrs = Vec::new();
        for (index,line) in file.lines().enumerate() {
            let line = line.unwrap();
            let split_line: Vec<&str> = line.split('\t').collect();
            match split_line[0].replace("chr", "").parse::<u8>() {
                Ok(result) =>{ 
                    let name      : u8        = result;
                    let length    : u32       = match split_line[1].parse() {
                        Ok(len) => len,
                        Err(e)  => return Err(ParseIntError(fai, name, e.to_string()))
                    };
                    debug!("Chromosome: {: <3} {: <10} {: <12}", index, name, length);
                    genome.add_chromosome(index, name, length);
                }
                Err(_) => skipped_chrs.push(split_line[0].to_string())
            };
        }
    
        if !skipped_chrs.is_empty(){
            warn!("\nSome chromosomes were skipped while parsing {}: {:#?}", fai, skipped_chrs);
        }
        Ok(genome)
    }

    pub fn add_chromosome(&mut self, index: usize, name: u8, length: u32) -> bool {
        self.0.insert(name, Chromosome::new(index, name, length)).is_none()
    }

    pub fn get_chr_mut(&mut self, name: &u8) -> Option<&mut Chromosome> {
        self.0.get_mut(name)
    }

    pub fn is_empty(&self) -> bool {
        self.0.values().all(|chr| chr.is_empty())
    }


    pub fn meiosis(&mut self, genetic_map: &GeneticMap) -> Gamete {
        let mut gamete = Gamete::default();
        for chr in self.0.values_mut() {
            if ! chr.is_empty(){
                let chromatid = chr.meiosis(genetic_map);
                gamete.add_chromatid(chr.name, chromatid);
            }
        }
        gamete
    }

}

/// Simply returns a default genome index in case the user did not provide a specific .fasta.fai file. 
impl Default for Genome {
    fn default() -> Self {
        warn!("No reference genome provided. Using default.");
        Genome(BTreeMap::from([
            (1, Chromosome::new( 0,  1,  249250621)),
            (2, Chromosome::new( 1,  2,  243199373)),
            (3, Chromosome::new( 2,  3,  198022430)),
            (4, Chromosome::new( 3,  4,  191154276)),
            (5, Chromosome::new( 4,  5,  180915260)),
            (6, Chromosome::new( 5,  6,  171115067)),
            (7, Chromosome::new( 6,  7,  159138663)),
            (8, Chromosome::new( 7,  8,  146364022)),
            (9, Chromosome::new( 8,  9,  141213431)),
            (10, Chromosome::new(9,  10, 135534747)),
            (11, Chromosome::new(10, 11, 135006516)),
            (12, Chromosome::new(11, 12, 133851895)),
            (13, Chromosome::new(12, 13, 115169878)),
            (14, Chromosome::new(13, 14, 107349540)),
            (15, Chromosome::new(14, 15, 102531392)),
            (16, Chromosome::new(15, 16,  90354753)),
            (17, Chromosome::new(16, 17,  81195210)),
            (18, Chromosome::new(17, 18,  78077248)),
            (19, Chromosome::new(18, 19,  59128983)),
            (20, Chromosome::new(19, 20,  63025520)),
            (21, Chromosome::new(20, 21,  48129895)),
            (22, Chromosome::new(21, 22,  51304566))
        ]))
    }
}