use std::borrow::Borrow;
use std::cmp::Ordering;

use crate::pileup::Nucleotide;
use genome::SNPCoord;

#[derive(Debug,Clone, Copy)]
pub struct Coordinate {
    pub chromosome: u8,
    pub position  : u32,
}

impl Coordinate {
    pub fn new(chromosome: u8, position: u32) -> Self {
        Self{chromosome, position}
    }
}

impl PartialEq<Coordinate> for Coordinate {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.position == other.position
    }
}

impl Eq for Coordinate {}

impl Ord for Coordinate {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.chromosome, self.position).cmp(&(other.chromosome, other.position))
    }
}

impl PartialOrd for Coordinate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


#[derive(Debug)]
pub struct Pwd {
    pub coordinate  : Coordinate,
    pub nucleotides : [Nucleotide; 2],
    pub pwd         : bool,
}

impl Pwd {
    pub fn new<'a>(coordinate: SNPCoord, random_nucl: &[&'a Nucleotide]) -> Self {
        Self {
            coordinate : Coordinate{chromosome: coordinate.chromosome, position: coordinate.position},
            nucleotides: [*random_nucl[0], *random_nucl[1]],
            pwd        : Self::check_pwd(&random_nucl),
        }
    }

    // Check if there is a pairwise difference.
    fn check_pwd<'a>(nuc: &[&'a Nucleotide]) -> bool {
        nuc[0].base != nuc[1].base
    }

    fn compute_avg_phred<'a>(&'a self) -> u8 {
        ( (self.nucleotides[0].phred + self.nucleotides[1].phred)/2 ) as u8
    }

    pub fn is_pwd(&self) -> bool {
        self.pwd
    }

    pub fn error_probs(&self) -> [f64; 2] {
        [
            self.nucleotides[0].error_prob(),
            self.nucleotides[1].error_prob()
        ]
    }
}

impl PartialEq<Pwd> for Pwd {
    fn eq(&self, other: &Self) -> bool { 
        self.coordinate == other.coordinate
    }
}

impl Eq for Pwd {}

impl Ord for Pwd {
    fn cmp(&self, other: &Self) -> Ordering {
        self.coordinate.cmp(&other.coordinate)
    }
}

impl PartialOrd for Pwd {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Borrow<Coordinate> for Pwd {
    fn borrow(&self) -> &Coordinate {
        &self.coordinate
    }
}