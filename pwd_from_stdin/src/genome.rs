use crate::jackknife::*;

use std::fs;
use std::hash::{Hash, Hasher};
use std::error::Error;
use std::io::{self, BufReader, BufRead};

use log::{info, warn};

#[derive(Debug)]
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
pub struct Chromosome {
    pub index  : usize,
    pub name   : u8,
    pub length : u32
}

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

pub fn fasta_index_reader(path: &String) -> Result<Vec<Chromosome>,Box<dyn Error>>  {
    info!("Parsing reference genome: {}", path);

    let mut genome: Vec<Chromosome> = Vec::new();
    let file = BufReader::new(fs::File::open(path).unwrap());

    let mut skipped_chrs = Vec::new();
    for (index,line) in file.lines().enumerate() {
        let line = line.unwrap();
        let split_line: Vec<&str> = line.split('\t').collect();
        match split_line[0].replace("chr", "").parse::<u8>() {
            Ok(result) =>{ 
                let name      : u8        = result;
                let length    : u32       = split_line[1].parse().unwrap();
                genome.push(Chromosome{index, name, length});
            }
            Err(_) => skipped_chrs.push(split_line[0].to_string())
        };
    }

    if !skipped_chrs.is_empty(){
        warn!("\nSome chromosomes were skipped while parsing {}: {:#?}", path, skipped_chrs);
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
