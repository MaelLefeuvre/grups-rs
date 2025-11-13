//use super::super::Individual;
pub use super::ComparisonError;

use fastrand::{self, Rng};
use located_error::LocatedOption;

use located_error::prelude::*;

use crate::pedigrees::{pedigree::individual::IndividualId};

// ------------------------------------------------------------------ //
// ----- Comparison
#[derive(Debug, Clone)]
pub struct PedComparison {
    pub label: String,
    pub pair: [IndividualId; 2],
    pwd: u32,
    overlap: u32,
}

impl PedComparison {
    pub fn new(label: &str, pair: [IndividualId; 2]) -> Self {
        Self {
            label: label.to_string(),
            pair,
            pwd: 0,
            overlap: 0,
        }
    }

    pub fn compare_alleles(&mut self, alleles: [[u8; 2]; 2], contam_rate: [f64; 2], contam_pop_af: [f64; 2], seq_error_rate: [f64; 2], rng: &mut Rng) -> Result<()> {
        use ComparisonError::CompareAllele;
        self.add_overlap();

        let random_sample0 = Self::simulate_observed_read(rng, contam_rate[0], contam_pop_af[0], seq_error_rate[0], alleles[0])
            .with_loc(||CompareAllele)?;
        let random_sample1 = Self::simulate_observed_read(rng, contam_rate[1], contam_pop_af[1], seq_error_rate[1], alleles[1])
            .with_loc(||CompareAllele)?;

        if random_sample0 != random_sample1 {
            self.add_pwd();
        }

        Ok(())
    }

    pub fn add_overlap(&mut self) {
        self.overlap += 1;
    }

    pub fn add_n_overlaps(&mut self, n: u32) {
        self.overlap += n;
    }

    pub fn add_pwd(&mut self) {
        self.pwd += 1;
    }

    /// Compute the average pairwise differences between the two individuals.
    /// 
    /// i.e. `self.pwd / self.overlap` 
    pub fn get_avg_pwd(&self) -> f64 {
        f64::from(self.pwd) / f64::from(self.overlap)
    }

    pub fn get_overlap(&self) -> u32 {
        self.overlap
    }

    pub fn get_sum_pwd(&self) -> u32 {
        self.pwd
    }

    #[inline]
    fn simulate_observed_read(rng: &mut fastrand::Rng, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Result<u8> {
        const SEQ_ERROR_CHOICES: [[u8; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        use ComparisonError::{SampleAllele, SimSeqError};
        // ---- Simulate modern human contamination. 
        let chosen_base: u8 = match rng.f64() < contam_rate {
            true  => match rng.f64() < contam_pop_af {
                true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                false => 0,  // otherwise, pick the reference allele.
            }
            false => *alleles.get(rng.usize(0..=1)).with_loc(||SampleAllele)?
        };

        // ---- Simulate sequencing error rate.
        // @ TODO: Find a smarter way to simulate error rates.
        if rng.f64() < seq_error_rate {
            let wrong_base: u8 = *SEQ_ERROR_CHOICES[chosen_base as usize].get(rng.usize(0..3)).with_loc(||SimSeqError)?;
            Ok(wrong_base)
        }
        else {
            Ok(chosen_base)
        }
    }
    
    /// Simulate `n` observed pileup reads from a set of alleles and given the provided contamination and sequencing parameters.
    /// # Arguments:
    /// - `n`             : Number of pileup observations to simulate.
    /// - `contam_rate`   : Modern human contamination rate required for the simulation.
    /// - `contam_pop_af` : allele frequency of the contaminating population for the current SNP coordinate.
    /// - `seq_error_rate`: sequencing error rate required for the simulation.
    /// - `alleles`       : size-two set of alleles of the pedigree individual for the current SNP coordinate.
    /// 
    /// # Note: 
    /// - this is a legacy function, which now wraps around the more performant [[`simulate_observed_read`]].
    /// - Keeping this for unit-testing purposes.
    #[cfg(test)]
    fn _simulate_observed_reads(n: u8, rng: &mut fastrand::Rng, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: [u8; 2]) -> Result<Vec<u8>> {
        // ---- Simulate n pileup observations.
        let mut reads = Vec::with_capacity(n as usize);
        for _ in 0..n {
            reads.push(Self::simulate_observed_read(rng, contam_rate, contam_pop_af, seq_error_rate, alleles)?);
        }
        Ok(reads)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pedigrees::pedigree::tests::common;
    use itertools::izip;

    fn get_expected_simulated_allele(haplo: u8, contam_rate: f64, contam_pop_af: f64) -> u8 {
        #![allow(clippy::cast_sign_loss,clippy::cast_possible_truncation)]
        (haplo + (contam_rate as u8 * (haplo + contam_pop_af as u8))) % 2
    }

    #[test]
    fn pwd_increment(){
        let mut pedigree = common::mock_pedcomparison();
        let comp = pedigree.comparisons.first_mut().expect("Comparison should be retrievable");
        assert_eq!(comp.pwd, 0);
        comp.add_pwd();
        assert_eq!(comp.pwd, 1);
    }

    #[test]
    fn overlap_increment(){
        let mut pedigree = common::mock_pedcomparison();
        let comp = pedigree.comparisons.first_mut().expect("Comparison should be retrievable");
        assert_eq!(comp.overlap, 0);
        comp.add_overlap();
        assert_eq!(comp.overlap, 1);
    }

    #[test]
    fn avg_pwd(){
        #![allow(clippy::float_cmp)]
        let n_iters = 10;
        let mut pedigree = common::mock_pedcomparison();
        let comp = pedigree.comparisons.first_mut().expect("Comparison should be retrievable");
        for pwd in 0..n_iters {
            comp.add_pwd();
            for overlap in 0..n_iters {
                comp.add_overlap();
                let want = f64::from(pwd + 1) / f64::from(overlap + (pwd*10) + 1);
                assert_eq!(comp.get_avg_pwd(), want);
            }
        }
    }


    #[test]
    fn simulate_observed_reads_contam() -> Result<()> {
        let binary_rates = [0.0, 1.0];
        let binary_alleles = [[0,0], [1,1]];
        let mut rng = fastrand::Rng::new();
        for contam_rate in binary_rates {
            for contam_pop_af in binary_rates {
                for alleles in binary_alleles {
                    let want = get_expected_simulated_allele(alleles[0], contam_rate, contam_pop_af);
                    let got = PedComparison::_simulate_observed_reads(1, &mut rng, contam_rate, contam_pop_af, 0.0, alleles)?;
                    assert_eq!(want, got[0]);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn allele_comparison() -> Result<()> {
        let mut rng = fastrand::Rng::new();

        let ref_alt = [0, 1];
        let binary_rates = [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]];
        for allele_ind_0 in ref_alt {
            for allele_ind_1 in ref_alt {
                for contam_rate in binary_rates {
                    for contam_pop_af in binary_rates {
                        let mut pedigree = common::mock_pedcomparison();
                        let comp = pedigree.comparisons.first().expect("Comparison should be retrievable");
                        let alleles = [[allele_ind_0, allele_ind_0], [allele_ind_1, allele_ind_1]];
                        comp.pair.into_iter().zip(alleles.iter()).for_each(|(ind, all)| {
                            pedigree.individuals.get_ind_mut(ind).expect("Individual should be retrievable").set_alleles(*all)
                        });

                        let comp = pedigree.comparisons.first_mut().expect("Comparison should be retrievable");
                        comp.compare_alleles(alleles, contam_rate, contam_pop_af, [0.0,0.0], &mut rng)?;

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
        Ok(())
    }
}