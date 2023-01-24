use super::{Contaminant, Pedigree};
use super::pedigree;

use std::{
    collections::HashMap,
    path::Path,
    ops::{Deref, DerefMut}
};

use grups_io::read::PanelReader;
use grups_io::read::SampleTag;

use located_error::prelude::*;

/// A vector of pedigree simulation replicates. This struct is generally assigned to a given Pileup-Comparison.
/// # Fields:
/// - `inner`       : vector of `Pedigree` simulation replicates.
/// - `contaminants`: Size-two array set of contaminating SampleTags, one vector for each compared pileup individual. 
///                   (see pedigree::Contaminant)
pub struct PedigreeReps{
    inner: Vec<Pedigree>,
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
    ///                         `sample_contam_tags[i] = contaminating individuals for pileup individual[i]
    /// -  `pair_indices`     : Size-two array, containing the pileup index of each pileup individual being compared.
    /// 
    /// # @ TODO:
    /// - `samples_contam_tags` should be an array. at the very least, this method should 
    ///    check if `samples_contam_tags.len()` == 2
    pub fn set_contaminants(&mut self, samples_contam_tags: &[Vec<SampleTag>], pair_indices: [usize; 2]) {
        let tags_0 = samples_contam_tags[pair_indices[0] % samples_contam_tags.len()].clone();
        let tags_1 = samples_contam_tags[pair_indices[1] % samples_contam_tags.len()].clone();
        self.contaminants = Some(Contaminant::new([tags_0, tags_1]))
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
    pub fn populate(&mut self, pedigree_path: &Path, pop: &String, panel: &PanelReader) -> Result<()> {
        let loc_msg = |ctxt: &str, i: usize| format!("While attempting to {ctxt} pedigree n°{i}");
        for i in 0..self.inner.capacity() {
            let mut pedigree = pedigree::parser::pedigree_parser(pedigree_path).with_loc(||loc_msg("parse", i))?;
            pedigree.set_tags(panel, pop, self.contaminants.as_ref()).with_loc(||loc_msg("set population tags of", i))?;
            pedigree.assign_offspring_strands().with_loc(||loc_msg("assign offspring strands of", i))?;
            self.inner.push(pedigree);
        }
        Ok(())
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
                let mut summary_statistics = sum_simulated_stats.get_mut(&comparison.label)
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
        let mut ordered_rels: Vec<_> = sum_simulated_stats.into_iter().map(|(rel, stats)| (rel, stats)).collect();
        ordered_rels.sort_by(|a, b| f64::total_cmp(&b.1.0, &a.1.0) );
        Ok(ordered_rels)
    }

    //// WIP: heterozygocity ratio
    //pub fn compute_average_het_ratio(&self) -> HashMap<String, f64>{
    //    let mut sum_simulated_stats = HashMap::new();
    //    for pedigree in self.iter() {
    //        // Sum the avg pwd of each replicate
    //        for comparison in pedigree.comparisons.iter() {
    //            *sum_simulated_stats.entry(comparison.label.to_owned()).or_insert(0.0) += comparison.get_heterozygocity_ratio()
    //        }
    //    }
    //    sum_simulated_stats
    //}

}

impl Deref for PedigreeReps {
    type Target = Vec<Pedigree>;
    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl DerefMut for PedigreeReps {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.inner
    }
}

impl std::fmt::Display for PedigreeReps {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.inner.iter().enumerate().fold(Ok(()), |_, (idx, pedigree)| {
            pedigree.comparisons.iter().fold(Ok(()), |result, comparison| {
                result.and_then(|_| writeln!(f, "{idx} - {}", comparison))
            })
        })
    }
}
