
use crate::jackknife::*;

use std::fs;
use std::hash::{Hash, Hasher};
use std::error::Error;
use std::io::{BufReader, BufRead};

use log::{info, warn, debug};

#[derive(Debug)]
/// A simple struct representing an SNPCoordinate position.
/// Note that :
///   - REF/ALT positions are optional, but are required if the user
///     requires known-sites filtration.
///   - SNPCoordinates implement Equality with JackknifeBlock structs.
///     See: `pwd_from_stdin::jackknife::JackknifeBlock`
///   - Hashable, but only in regards to chr and pos.
///     -> don't insert alternate variations.
pub struct SNPCoord {
    pub chromosome : u8,
    pub position   : u32,
    pub reference  : Option<char>,
    pub alternate  : Option<char>,
}


impl PartialEq<SNPCoord> for SNPCoord {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.position == other.position
    }
}

impl PartialEq<JackknifeBlock> for SNPCoord {
    fn eq(&self, other: &JackknifeBlock) -> bool {
        other.chromosome == self.chromosome && self.position >= other.range.start && self.position < other.range.end
    }
}

impl Eq for SNPCoord {}

impl Hash for SNPCoord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.position.hash(state);
    }
}


#[derive(Debug)]
/// A simple struct representing a chromosome index. This is mainly used to compute Jackknife Blocks
pub struct Chromosome {
    pub index  : usize,
    pub name   : u8,
    pub length : u32
}

/// Simply returns a default genome index in case the user did not provide a specific .fasta.fai file. 
pub fn default_genome() -> Vec<Chromosome> {
    warn!("No reference genome provided. Using default.");
    vec![
        Chromosome{index:  0, name:  1, length: 249250621},
        Chromosome{index:  1, name:  2, length: 243199373},
        Chromosome{index:  2, name:  3, length: 198022430},
        Chromosome{index:  3, name:  4, length: 191154276},
        Chromosome{index:  4, name:  5, length: 180915260},
        Chromosome{index:  5, name:  6, length: 171115067},
        Chromosome{index:  6, name:  7, length: 159138663},
        Chromosome{index:  7, name:  8, length: 146364022},
        Chromosome{index:  8, name:  9, length: 141213431},
        Chromosome{index:  9, name: 10, length: 135534747},
        Chromosome{index: 10, name: 11, length: 135006516},
        Chromosome{index: 11, name: 12, length: 133851895},
        Chromosome{index: 12, name: 13, length: 115169878},
        Chromosome{index: 13, name: 14, length: 107349540},
        Chromosome{index: 14, name: 15, length: 102531392},
        Chromosome{index: 15, name: 16, length:  90354753},
        Chromosome{index: 16, name: 17, length:  81195210},
        Chromosome{index: 17, name: 18, length:  78077248},
        Chromosome{index: 18, name: 19, length:  59128983},
        Chromosome{index: 19, name: 20, length:  63025520},
        Chromosome{index: 20, name: 21, length:  48129895},
        Chromosome{index: 21, name: 22, length:  51304566}
    ]
}




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
            Self::InvalidFileFormatError(ext)  => write!(f, "Invalid fasta format: expected [.fasta|.fa], got {}", ext),
        }
    }
}

use std::path::Path;
/// Read a `.fasta.fai` file and parse it into a vector of Chromosome structs.
/// 
/// TODO: - convert 'path' variable to &str
///       - Check for .fai.fai extension doubles.
pub fn fasta_index_reader(path: &str) -> Result<Vec<Chromosome>,FastaIndexReaderError>  {
    info!("Parsing reference genome: {}", path);

    let fai = match Path::new(path).extension().unwrap().to_str() { // Append '.fai' if it is not there 
        Some("fai")          => path.clone().to_string(),                          
        Some("fasta" | "fa") => format!("{}{}", path.clone(), ".fai"),
        Some(other)          => return Err(FastaIndexReaderError::InvalidFileFormatError(other.to_string())),
        None                 => return Err(FastaIndexReaderError::InvalidFileFormatError("None".to_string())),
        
    };
    debug!("fasta.fai file : {}", fai);


    let mut genome: Vec<Chromosome> = Vec::new();
    let file = BufReader::new(match fs::File::open(fai.clone()) {
        Ok(file) => file,
        Err(err) => return Err(FastaIndexReaderError::FileNotFoundError(fai, err.to_string())),
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
                    Err(e)  => return Err(FastaIndexReaderError::ParseIntError(fai, name, e.to_string()))
                };
                debug!("Chromosome: {: <3} {: <10} {: <12}", index, name, length);
                genome.push(Chromosome{index, name, length});
            }
            Err(_) => skipped_chrs.push(split_line[0].to_string())
        };
    }

    if !skipped_chrs.is_empty(){
        warn!("\nSome chromosomes were skipped while parsing {}: {:#?}", fai, skipped_chrs);
    }
    Ok(genome)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_full_equality() {
        let chromosome1 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        let chromosome2 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }

    #[test]
    fn snp_partial_equality() {
        let chromosome1 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('T'), alternate: Some('G')};
        let chromosome2 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }
}
