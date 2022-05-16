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

use genome::{
    Genome,
};


pub struct PedigreeReps{
    inner: Vec<Pedigree>,
    pub contaminants: Option<Contaminant> 
}

impl PedigreeReps {
    pub fn with_capacity(n: usize) -> Self {
        Self{inner: Vec::with_capacity(n), contaminants: None}
    }
    
    pub fn set_contaminants(&mut self, samples_contam_tags: &Vec<Vec<SampleTag>>, pair_indices: [usize; 2]) {
        let tags_0 = samples_contam_tags[pair_indices[0] % samples_contam_tags.len()].clone(); // Wrap if contam_set.len() < comparison.len()
        let tags_1 = samples_contam_tags[pair_indices[1] % samples_contam_tags.len()].clone();
        self.contaminants = Some(Contaminant::new([tags_0, tags_1]))
    }

    pub fn populate(&mut self, pedigree_path: &Path, pop: &String, panel: &VCFPanelReader, genome: &Genome) -> Result<(), Box<dyn Error>> {
        for _ in 0..self.inner.capacity() {
            let mut pedigree = pedigree::io::pedigree_parser(pedigree_path, genome).unwrap();
            pedigree.set_tags(panel, pop, self.contaminants.as_ref());
            pedigree.assign_offspring_strands()?;
            self.inner.push(pedigree);
        }
        Ok(())
    }

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
