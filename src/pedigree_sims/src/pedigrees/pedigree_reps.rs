use super::pedigree::parser::PedigreeBuilder;
use super::{Contaminant, Pedigree};
use crate::pedigrees::constants::REPLICATE_ID_FORMAT_LEN;

use std::{
    collections::HashMap,
    path::Path,
    ops::{Deref, DerefMut},
    fmt::{self, Formatter, Display}
};

use grups_io::read::PanelReader;
use grups_io::read::SampleTag;

use located_error::prelude::*;

/// A vector of pedigree simulation replicates. This struct is generally assigned to a given Pileup-Comparison.
/// # Fields:
/// - `inner`       : vector of `Pedigree` simulation replicates.
/// - `contaminants`: Size-two array set of contaminating SampleTags, one vector for each compared pileup individual. (see pedigree::Contaminant)
pub struct PedigreeReps{
    pub inner: Vec<Pedigree>,
    pub contaminants: Option<Contaminant> 
}

impl PedigreeReps {
    /// Instantiate an empty vector of pedigrees with a preset capacity.
    /// # Arguments:
    /// - `n`: allocated size of the `self.inner` vector.
    pub fn with_capacity(n: usize) -> Self {
        Self{inner: Vec::with_capacity(n), contaminants: None}
    }
    
    /// Instantiate and assign a new `Contaminant` object to `self.contaminants`
    /// # Arguments
    /// - `sample_contam_tags`: Size-two vector of vectors of contaminating SampleTags.
    ///   - `sample_contam_tags[i] = contaminating individuals for pileup individual[i]
    /// -  `pair_indices`     : Size-two array, containing the pileup index of each pileup individual being compared.
    /// 
    /// # @ TODO:
    /// - `samples_contam_tags` should be an array. at the very least, this method should 
    ///   check if `samples_contam_tags.len()` == 2
    pub fn set_contaminants(&mut self, samples_contam_tags: &[Vec<SampleTag>], pair_indices: [usize; 2]) {
        let tags_0 = samples_contam_tags[pair_indices[0] % samples_contam_tags.len()].clone();
        let tags_1 = samples_contam_tags[pair_indices[1] % samples_contam_tags.len()].clone();
        self.contaminants = Some(Contaminant::new([tags_0, tags_1]));
    }

    /// Instantiate and insert `n` pedigrees within `self.inner` (`n`, being the capacity of `self.inner)
    /// # Arguments
    /// - `pedigree_path`: path leading to the pedigree definition file. 
    /// - `pop`          : name of the pedigree's desired (super-)population of origin
    /// - `panel`        : input samples definition file, in the form of a `PanelReader`.
    ///
    /// # Errors:
    /// - if any of the offspring `self.parents` field is set to None
    /// - if pop label is invalid.
    /// - returns `std::io::result::InvalidData` if an error occurs while parsing the pedigree definition file.
    pub fn populate(&mut self, pedigree_path: &Path) -> Result<()> {
        let loc_msg = |ctxt: &str, i: usize| format!("While attempting to {ctxt} pedigree n°{i}");
        let pedigree_builder = PedigreeBuilder::new(pedigree_path)
            .with_loc(|| format!("While attempting to read the pedigree definition file {}", pedigree_path.display()))?;
        for i in 0..self.inner.capacity() {
            let pedigree = pedigree_builder.build().with_loc(||loc_msg("parse", i))?;
            self.inner.push(pedigree);
        }
        Ok(())
    }

    pub fn assign_offspring_strands(&mut self) -> Result<()> {
        self.inner.iter_mut().enumerate().try_for_each(|(i, pedigree)|{
            pedigree.assign_offspring_strands().with_loc(||format!("While attempting to assign offspring strands of pedigree n°{i}"))
        })
    }
    pub fn set_founder_tags(&mut self, panel: &PanelReader, pop: &String) -> Result<()>{
        self.inner.iter_mut().enumerate().try_for_each(|(i, pedigree)| {
            pedigree.set_founder_tags(panel, pop, self.contaminants.as_ref()).with_loc(||format!("While attempting to set population tags of pedigree n°{i}"))
        })
    }

    pub fn assign_random_sex(&mut self) -> Result<()> {
        self.inner.iter_mut().enumerate().try_for_each(|(i, pedigree)| {
            pedigree.assign_random_sex().with_loc(|| format!("While attempting to randomly assign sex of pedigree n°{i}"))
        })
    }

    /// Aggregate the sum of all simulated pairwise differences for each pedigree comparison.
    /// Returns a hashmap with Key = comparison_label | Value = ( avg_avg_pwd, (avg_avg_pwd) )
    /// # @TODO: This data structure is error prone and should be converted to a struct or named tuple.
    pub fn compute_sum_simulated_stats(&self) -> Result<Vec<(String, (f64,f64))>> {
        let mut sum_simulated_stats = HashMap::new();
        for pedigree in self.iter() {
            // Sum the avg pwd of each replicate.
            for comparison in pedigree.comparisons.iter() {
                sum_simulated_stats.entry(comparison.label.to_owned()).or_insert((0.0, 0.0)).0 += comparison.get_avg_pwd()
            }
        }

        // ---- Compute the sum-squared for sample variance for each comparison.
        for pedigree in self.iter() {
            for comparison in pedigree.comparisons.iter() {
                // ---- Access summary statistics 
                let summary_statistics = sum_simulated_stats.get_mut(&comparison.label)
                    .with_loc(|| format!("While computing simulation summary statistics:\
                        Attempting to access a missing summary statistic using the comparison label {}",
                        comparison.label
                    ))?;

                // ---- Compute the avg + sum of squares for variance / std err estimation.
                let simulation_avg = summary_statistics.0 / self.len() as f64;
                summary_statistics.1 += (comparison.get_avg_pwd() - simulation_avg).powf(2.0);
            }
        }

        // Divide by the length of the inner vector to get the avg_avg, and standard deviation.
        let divisor = self.inner.len() as f64;
        for stat in sum_simulated_stats.values_mut() {
            stat.0 /= divisor;                          // avg of averages.
            stat.1 = (stat.1 / (divisor - 1.0)).sqrt(); //std.dev.
        }

        // Sort according to decreasing values of Avg(pwd)
        let mut ordered_rels: Vec<_> = sum_simulated_stats.into_iter().collect();
        ordered_rels.sort_by(|a, b| f64::total_cmp(&b.1.0, &a.1.0) );
        Ok(ordered_rels)
    }


    /// Apply `add_n_overlap` to every comparison of every contained pedigree.
    #[inline]
    pub fn add_non_informative_snps(&mut self, n: u32) {
        for pedigree in self.inner.iter_mut() {
            pedigree.comparisons.iter_mut().for_each(|comp| comp.add_n_overlaps(n));
        }
    }

    pub fn all_sex_assigned(&self) -> bool {
        self.inner.iter().all(|ped| ped.all_sex_assigned())
    }
}

impl Deref for PedigreeReps {
    type Target = [Pedigree];
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl DerefMut for PedigreeReps {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl Display for PedigreeReps {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        self.inner.iter().enumerate().try_fold((), |_, (idx, pedigree)| {
            pedigree.comparisons.iter().try_fold((), |_, comparison| {
                writeln!(f, "{idx: <REPLICATE_ID_FORMAT_LEN$} - {comparison}")
            })
        })
    }
}