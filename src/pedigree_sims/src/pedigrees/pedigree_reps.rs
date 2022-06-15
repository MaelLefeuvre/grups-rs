use super::{Contaminant, Pedigree};
use super::pedigree;

use std::{
    collections::HashMap,
    path::Path,
    error::Error,
    ops::{Deref, DerefMut}
};

use crate::io::vcf::{
    SampleTag,
    reader::VCFPanelReader,
};

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
    pub fn set_contaminants(&mut self, samples_contam_tags: &Vec<Vec<SampleTag>>, pair_indices: [usize; 2]) {
        let tags_0 = samples_contam_tags[pair_indices[0] % samples_contam_tags.len()].clone();
        let tags_1 = samples_contam_tags[pair_indices[1] % samples_contam_tags.len()].clone();
        self.contaminants = Some(Contaminant::new([tags_0, tags_1]))
    }

    /// Instantiate and insert `n` pedigrees within `self.inner` (`n`, being the capacity of `self.inner)
    /// # Arguments
    /// - `pedigree_path`: path leading to the pedigree definition file. 
    /// - `pop`          : name of the pedigree's desired (super-)population of origin
    /// - `panel`        : input samples definition file, in the form of a `VCFPanelReader`.
    ///
    /// # Errors:
    /// - if any of the offspring `self.parents` field is set to None
    /// - if pop label is invalid.
    /// - returns `std::io::result::InvalidData` if an error occurs while parsing the pedigree definition file.
    pub fn populate(&mut self, pedigree_path: &Path, pop: &String, panel: &VCFPanelReader) -> Result<(), Box<dyn Error>> {
        for _ in 0..self.inner.capacity() {
            let mut pedigree = pedigree::io::pedigree_parser(pedigree_path)?;
            pedigree.set_tags(panel, pop, self.contaminants.as_ref())?;
            pedigree.assign_offspring_strands()?;
            self.inner.push(pedigree);
        }
        Ok(())
    }

    /// Aggregate the sum of all simulated pairwise differences for each pedigree comparison.
    /// Returns a hashmap with Key = comparison_label | Value = sum( avg_pwd )
    pub fn compute_sum_simulated_pwds(&self) -> HashMap<String, f64> {
        let mut sum_simulated_pwds = HashMap::new();
        for pedigree in self.iter() {
            for comparison in pedigree.comparisons.iter() {
                *sum_simulated_pwds.entry(comparison.label.to_owned()).or_insert(0.0) += comparison.get_avg_pwd()
            }
        }
        sum_simulated_pwds
    }
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
                result.and_then(|_| writeln!(f, "{idx} {}", comparison))
            })
        })
    }
}
