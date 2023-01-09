use std::borrow::Borrow;
use std::cmp::Ordering;

use crate::pileup::{Nucleotide, Line};
use genome::SNPCoord;

use super::Individual;

use genome::coordinate::Coordinate;

#[derive(Debug)]
pub struct Pwd {
    pub coordinate  : Coordinate,
    pub phreds      : [f64; 2],
    pub pwd         : f64,
    observations    : u32,
}

impl Pwd {
    #[must_use]
    pub fn initialize(coordinate: SNPCoord) -> Self {
        Self{
            coordinate  : Coordinate{chromosome: coordinate.coordinate.chromosome, position: coordinate.coordinate.position},
            phreds      : [0.0,0.0],
            pwd         : 0.0,
            observations: 0
        }
    }
    
    #[must_use]
    pub fn one(coordinate: SNPCoord, random_nucl: &[&Nucleotide]) -> Self {
        Self {
            coordinate  : Coordinate{chromosome: coordinate.coordinate.chromosome, position: coordinate.coordinate.position},
            //nucleotides: [*random_nucl[0], *random_nucl[1]],
            phreds      : [
                f64::from(random_nucl[0].phred),
                f64::from(random_nucl[1].phred)
            ],
            pwd         : Self::check_pwd(random_nucl),
            observations: 1,
        }
    }

    #[must_use]
    pub fn deterministic_self(line: &Line, pair: &[Individual; 2]) -> Self {
        use itertools::Itertools;
        let (mut pwd, mut counter) = (0.0, 0.0);
        let mut phreds = [0.0, 0.0];
        // let mut hom_alt_sum = 0.0; // WIP: heterozygocity ratio
        for nucs in line.individuals[pair[0].index].nucleotides.iter().combinations(2) {
            if nucs[0].base != nucs[1].base {
                pwd += 1.0;
            } //else if nucs[0].base != '.' { // WIP: heterozygocity ratio
            //    hom_alt_sum += 1.0
            //}
            phreds[0] += f64::from(nucs[0].phred);
            phreds[1] += f64::from(nucs[1].phred);
            counter += 1.0; 
        }
        
        let coordinate = Coordinate{chromosome: line.coordinate.chromosome, position: line.coordinate.position};
        let phreds = [phreds[0]/counter , phreds[1]/counter];
        let pwd = pwd/counter ;

        Self { coordinate, phreds, pwd, observations: 1 }
    }

    #[must_use]
    pub fn deterministic_pairwise(line: &Line, pair: &[Individual; 2]) -> Self {
        // Breaks if self.comparison == true
        let set0 = line.individuals[pair[0].index].observation_set();
        let set1 = line.individuals[pair[1].index].observation_set();
        let phreds = [set0.1, set1.1];
        let mut prob_pwd = 0.0;

        //let mut hom_alt_sum = 0.0;

        for (base0, prob0) in &set0.0 {
            for (base1, prob1) in &set1.0 {
                if base0 != base1 {
                    prob_pwd += prob0 * prob1;
                    
                } //else if *base0 != '.' { // WIP: heterozygocity ratio
                //    hom_alt_sum += 1.0
                //}
            }
        }

        let coordinate = Coordinate{chromosome: line.coordinate.chromosome, position: line.coordinate.position};
        Self { coordinate, phreds, pwd: prob_pwd, observations: 1 }
    }

    pub fn update(&mut self, random_nucl: &[&Nucleotide]) {
        self.pwd += Self::check_pwd(random_nucl);
        self.update_phreds(random_nucl);
        self.observations += 1;
    }

    fn update_phreds(&mut self, random_nucl: &[&Nucleotide]) {
        self.phreds[0] += f64::from(random_nucl[0].phred);
        self.phreds[1] += f64::from(random_nucl[1].phred);
    }
    
    // Check if there is a pairwise difference.
    fn check_pwd(nuc: &[&Nucleotide]) -> f64 {
        f64::from(u8::from(nuc[0].base != nuc[1].base))
    }

    #[must_use]
    pub fn avg_local_pwd(&self) -> f64 {
        self.pwd / f64::from(self.observations)
    }

    #[must_use]
    pub fn compute_avg_phred(&self) -> f64 {
        (self.phreds[0] + self.phreds[1]) / 2.0 / f64::from(self.observations)
    }

    #[must_use]
    pub fn error_probs(&self) -> [f64; 2] {
        [
            f64::powf(10.0, -1.0 * (self.phreds[0]) / 10.0),
            f64::powf(10.0, -1.0 * (self.phreds[1]) / 10.0)
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

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;
    use crate::pileup::Line;

    #[test]
    pub fn deterministic_pairwise() -> Result<(), Box<dyn Error>> {
        let raw_line="22\t51057923\tC\t6\tTTTTtt\tJEJJEE\t4\t.,TT\tJJJJ";
        let line = Line::new(raw_line, true)?;
        let pair = [
            Individual::new(None, 0, 1),
            Individual::new(None, 1, 1),
        ];
        let pwd = Pwd::deterministic_pairwise(&line, &pair);
        assert_eq!(pwd.pwd, 0.5);
        assert_eq!(pwd.phreds, [38.5, 41.0]);
        Ok(())
    }

    #[test]
    pub fn deterministic_self() -> Result<(), Box<dyn Error>> {
        let raw_line="22\t51057923\tC\t4\tTTT...\tJJJEEE";
        let line = Line::new(raw_line, true)?;
        let pair = [
            Individual::new(None, 0, 2),
            Individual::new(None, 0, 2),
        ];
        let pwd = Pwd::deterministic_self(&line, &pair);
        assert_eq!(pwd.pwd, 0.6);
        assert_eq!(pwd.phreds, [40.0, 37.0]);
        Ok(())
    }

}