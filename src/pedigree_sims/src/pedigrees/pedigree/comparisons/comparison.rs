use std::{
    error::Error,
    rc::Rc,
    cell::RefCell
};

use super::super::Individual;

use crate::io::vcf::SampleTag; 



use rand::Rng;
use rand::seq::SliceRandom;

#[derive(Debug, Clone)]
pub struct PedComparison {
    pub label            : String,
    pair             : (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>),
    pwd              : u32,
    overlap          : u32,
    _self_comparison : bool,
}

impl PedComparison {
    pub fn new(label: String, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>), self_comparison: bool) -> PedComparison {
        let pair = Self::format_pair(pair);
        PedComparison{label, pair, _self_comparison: self_comparison, pwd: 0, overlap: 0 }
    }

    pub fn add_pwd(&mut self){
        self.pwd += 1;
    }

    pub fn add_overlap(&mut self){
        self.overlap +=1;
    }

    pub fn get_avg_pwd(&self) -> f64 {
        self.pwd as f64 / self.overlap as f64
    }
    
    fn format_pair(pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>) {
        (Rc::clone(pair.0), Rc::clone(pair.1))
    }

    pub fn compare_alleles(&mut self, contam_rate: [f64; 2], contam_pop_af: [f64; 2], seq_error_rate: [f64; 2]) -> Result<(), Box<dyn Error>> {
        self.add_overlap();
        let random_sample0 = Self::simulate_observed_reads(1, contam_rate[0], contam_pop_af[0], seq_error_rate[0], self.pair.0.borrow().get_alleles()?);
        let random_sample1 = Self::simulate_observed_reads(1, contam_rate[1], contam_pop_af[1], seq_error_rate[1], self.pair.1.borrow().get_alleles()?);
        if random_sample0 != random_sample1 {
            self.add_pwd();
        }
        Ok(())
    }

    fn simulate_observed_reads(n: u8, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Vec<u8> {
        let mut reads = Vec::with_capacity(n as usize);
        let mut rng = rand::thread_rng();

        // Simulate contamination.
        for _ in 0..n {
            let chosen_base: u8 = match rng.gen::<f64>() < contam_rate {
                true  => match rng.gen::<f64>() < contam_pop_af {
                    true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                    false => 0,  // otherwise, pick the reference allele.
                }
                false => *alleles.choose(&mut rng).unwrap(),
            };

            // Simulate sequencing error rate.
            let seqerror_choices: Vec<[u8; 3]> = vec![[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
            if rng.gen::<f64>() < seq_error_rate {
                let wrong_base: u8 = *seqerror_choices[chosen_base as usize].choose(&mut rng).unwrap();
                reads.push(wrong_base);
            }
            else {
                reads.push(chosen_base);
            }
        }
        reads
    }
}

impl std::fmt::Display for PedComparison {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let default_tag=SampleTag::new("None", None);
        write!(f, " {: <20} {: <8} {: <8} {: <8} {: <8} {: >9} {: >9} {: <12.6}",
            self.label,
            self.pair.0.borrow().label,
            self.pair.1.borrow().label,
            self.pair.0.borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pair.1.borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pwd,
            self.overlap,
            self.get_avg_pwd()
        )
    }
}
