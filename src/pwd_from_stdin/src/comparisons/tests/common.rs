
use std::error::Error;

use super::super::comparison::Comparison;
use crate::comparisons::Individual;
use crate::pileup::Pileup;
use genome::Genome;
use genome::snp::Allele;
use rand::prelude::SliceRandom;

pub const MOCK_IND_INDICES: [usize; 2] = [0, 1];
pub const MOCK_IND_MINDEPTH: [u16; 2] = [2; 2];
pub const MOCK_BLOCKSIZE: u32 = 1000;


pub fn mock_comparison(self_comparison: bool) -> Comparison {
    let mockind_0 = Individual::new(None, MOCK_IND_INDICES[0], MOCK_IND_MINDEPTH[0]);
    let mockind_1 = match self_comparison {
        false  => Individual::new(None, MOCK_IND_INDICES[1], MOCK_IND_MINDEPTH[1]),
        true   => mockind_0.clone(),
    };

    let default_genome = Genome::default();
    Comparison::new([mockind_0, mockind_1], self_comparison, &default_genome, MOCK_BLOCKSIZE)
}


fn mock_pileup_strings(depth: u16, min_qual: u8) -> Result<(String, String), Box<dyn Error>> {
    let mut rng = rand::thread_rng();

    let nucleotides = vec!['A', 'C', 'G', 'T'];
    let quals: Vec<char> = std::str::from_utf8( &(min_qual+33..33+33).collect::<Vec<u8>>())?.chars().collect();

    let mut bases = String::new();
    let mut scores = String::new();
    for _ in 0..depth+1 {
        bases.push(*nucleotides.choose(&mut rng).ok_or("Empty nucleotide vector in 'mock_pileup_strings'")?);
        scores.push(*quals.choose(&mut rng).ok_or("Empty quals vector in 'mock_pileup_strings'")?);
    }

    Ok((bases, scores))
}


pub fn mock_pileups(min_depths: &[u16], min_qual: u8, ignore_dels: bool) -> Result<Vec<Pileup>, Box<dyn Error>> {
    let mut pileups = Vec::new();
    for depth in min_depths {
        let (bases, scores) = mock_pileup_strings(*depth, min_qual)?;
        let pileup = Pileup::new(Allele::N, *depth, bases.as_str(), scores.as_str(), ignore_dels)?;
        pileups.push(pileup);
    }

    Ok(pileups)
}