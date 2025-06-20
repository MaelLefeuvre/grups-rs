use crate::pileup::Line;
use genome::Nucleotide;

use super::Individual;

use genome::coordinate::{Coordinate, derive::{Coord, CoordBorrow, CoordEq, CoordHash, CoordOrd}};

#[derive(Debug, Coord, CoordEq, CoordOrd, CoordHash, CoordBorrow)]
pub struct Pwd {
    pub coordinate  : Coordinate,
    pub phred_sums  : [f64; 2],
    pub pwd         : f64,
    observations    : u32,
}

impl Pwd {
    #[must_use]
    pub fn initialize(coordinate: Coordinate) -> Self {
        Self{
            coordinate,
            phred_sums  : [0.0,0.0],
            pwd         : 0.0,
            observations: 0
        }
    }
    
    #[must_use]
    pub fn one(coordinate: Coordinate, random_nucl: &[&Nucleotide]) -> Self {
        Self {
            coordinate,
            phred_sums  : [0, 1].map(|i| f64::from(random_nucl[i].phred.score()) ),
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
            [0, 1].into_iter().for_each(|i| phreds[i] += f64::from(nucs[i].phred.score()));
            counter += 1.0; 
        }
        
        let phred_sums = [phreds[0]/counter , phreds[1]/counter];
        let pwd = pwd/counter ;

        Self { coordinate: line.coordinate, phred_sums, pwd, observations: 1 }
    }

    #[must_use]
    pub fn deterministic_pairwise(line: &Line, pair: &[Individual; 2]) -> Self {
        // Breaks if self.comparison == true
        let observation_sets = [0, 1].map(|i| line.individuals[pair[i].index].observation_set());
        
        let mut prob_pwd = 0.0;

        //let mut hom_alt_sum = 0.0;
        for (base0, prob0) in &observation_sets[0].0 {
            for (base1, prob1) in &observation_sets[1].0 {
                if base0 != base1 {
                    prob_pwd += prob0 * prob1;
                    
                } //else if *base0 != '.' { // WIP: heterozygocity ratio
                //    hom_alt_sum += 1.0
                //}
            }
        }
        let coordinate = Coordinate{chromosome: line.coordinate.chromosome, position: line.coordinate.position};
        let phred_sums = observation_sets.map(|set| set.1);
        Self { coordinate, phred_sums, pwd: prob_pwd, observations: 1 }
    }

    pub fn update(&mut self, random_nucl: &[&Nucleotide]) {
        self.pwd += Self::check_pwd(random_nucl);
        self.update_phreds(random_nucl);
        self.observations += 1;
    }

    fn update_phreds(&mut self, random_nucl: &[&Nucleotide]) {
        [0,1].into_iter().for_each(|i| self.phred_sums[i] += f64::from(random_nucl[i].phred.score()) );
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
        f64::midpoint(self.phred_sums[0], self.phred_sums[1]) / f64::from(self.observations)
    }

    #[must_use]
    pub fn error_probs(&self) -> [f64; 2] {
        // @ TODO: this is very wrong, since we're re-implementing a phred method.
        [0, 1].map(|i|  f64::powf(10.0, -1.0 * (self.phred_sums[i]) / 10.0) )
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]
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
        assert_eq!(pwd.phred_sums, [38.5, 41.0]);
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
        assert_eq!(pwd.phred_sums, [40.0, 37.0]);
        Ok(())
    }

}