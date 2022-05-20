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
    pub label        : String,
    pair             : [Rc<RefCell<Individual>>; 2],
    pwd              : u32,
    overlap          : u32,
    _self_comparison : bool,
}

impl PedComparison {
    pub fn new(label: &str, pair: [&Rc<RefCell<Individual>>; 2], self_comparison: bool) -> PedComparison {
        let pair = Self::format_pair(pair);
        PedComparison{label: label.to_string(), pair, _self_comparison: self_comparison, pwd: 0, overlap: 0 }
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
    
    fn format_pair(pair: [&Rc<RefCell<Individual>>; 2]) ->  [Rc<RefCell<Individual>>; 2] {
        [Rc::clone(pair[0]), Rc::clone(pair[1])]
    }

    pub fn compare_alleles(&mut self, contam_rate: [f64; 2], contam_pop_af: [f64; 2], seq_error_rate: [f64; 2]) -> Result<(), Box<dyn Error>> {
        self.add_overlap();
        let random_sample0 = Self::simulate_observed_reads(1, contam_rate[0], contam_pop_af[0], seq_error_rate[0], self.pair[0].borrow().get_alleles()?);
        let random_sample1 = Self::simulate_observed_reads(1, contam_rate[1], contam_pop_af[1], seq_error_rate[1], self.pair[1].borrow().get_alleles()?);
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
        write!(f, "{: <20} - {: <8} - {: <8} - {: <8} - {: <8} - {: >9} - {: >9} - {: <12.6}",
            self.label,
            self.pair[0].borrow().label,
            self.pair[1].borrow().label,
            self.pair[0].borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pair[1].borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pwd,
            self.overlap,
            self.get_avg_pwd()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pedigrees::pedigree::tests::common;
    use itertools::izip;

    fn get_expected_simulated_allele(haplo: u8, contam_rate: f64, contam_pop_af: f64) -> u8 {
        (haplo + (contam_rate as u8 * (haplo + contam_pop_af as u8))) % 2
    }

    #[test]
    fn pwd_increment(){
        let mut comp = common::mock_pedcomparison();
        assert_eq!(comp.pwd, 0);
        comp.add_pwd();
        assert_eq!(comp.pwd, 1);
    }

    #[test]
    fn overlap_increment(){
        let mut comp = common::mock_pedcomparison();
        assert_eq!(comp.overlap, 0);
        comp.add_overlap();
        assert_eq!(comp.overlap, 1);
    }

    #[test]
    fn avg_pwd(){
        let n_iters = 10;
        let mut comp = common::mock_pedcomparison();
        for pwd in 0..n_iters {
            comp.add_pwd();
            for overlap in 0..n_iters {
                comp.add_overlap();
                let want = (pwd+1) as f64 / (overlap + (pwd*10) +1) as f64;
                assert_eq!(comp.get_avg_pwd(), want);
            }
        }
    }


    #[test]
    fn simulate_observed_reads_contam() {
        let binary_rates = [0.0, 1.0];
        let binary_alleles = [[0,0], [1,1]];
        for contam_rate in binary_rates {
            for contam_pop_af in binary_rates {
                for alleles in binary_alleles {
                    let want = get_expected_simulated_allele(alleles[0], contam_rate, contam_pop_af);
                    let got = PedComparison::simulate_observed_reads(1, contam_rate, contam_pop_af, 0.0, alleles);
                    assert_eq!(want, got[0]);
                }
            }
        }
    }

    #[test]
    fn allele_comparison(){
        let ref_alt = [0, 1];

        let binary_rates = [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]];

        for allele_ind_0 in ref_alt {
            for allele_ind_1 in ref_alt {
                for contam_rate in binary_rates {
                    for contam_pop_af in binary_rates {
                        let mut comp = common::mock_pedcomparison();
                        let alleles = [[allele_ind_0, allele_ind_0], [allele_ind_1, allele_ind_1]];
                        comp.pair.iter().zip(alleles.iter()).for_each(|(ind, all)| ind.borrow_mut().set_alleles(*all));
                        comp.compare_alleles(contam_rate, contam_pop_af, [0.0,0.0]).unwrap();

                        let mut want = [0, 0];
                        izip!(&mut want, [allele_ind_0, allele_ind_1], contam_rate, contam_pop_af)
                            .for_each( |(res, all, c_rate, c_af)| {
                                *res += get_expected_simulated_allele(all, c_rate, c_af);
                            });

                        assert_eq!(want[0] != want[1], comp.get_avg_pwd() != 0.0);
                    }
                }
            }
        }
    }

    #[test]
    fn display() {
        let comp = common::mock_pedcomparison();
        println!("{comp}")
    }


}