use std::collections::BTreeMap;

use located_error::prelude::*;
use grups_io::read::{
    PanelReader,
    genotype_reader::GenotypeReader
};

use genome::Sex;

use fastrand;

use std::{
    cell::{RefCell, RefMut},
    rc::Rc,
};

#[cfg(test)] mod tests ; 

mod individual;
use individual::Individual;

pub mod parser;

mod contaminant;
pub use contaminant::Contaminant;

mod comparisons;
use comparisons::PedComparisons;
use comparisons::PedComparison;


pub mod pedparam;
use pedparam::PedigreeParams;

mod error;
pub use error::PedigreeError;

/// A Pedigree simulation replicate.
/// # Fields:
/// - `individuals`: BTreeMap containing all members of the pedigree (founders and offspring)
///                  - Key  : (String)                   - Label of the individual
///                  - Value: (Rc<RefCell<Individual>>>) - Reference to an individual, with interior mutability.
/// - `comparisons`: `PedComparison` struct, containing all the user-requested kinship scenarios.
/// - `params`     : `PedigreeParams` wrapper struct, containing the set of parameters (contam_rate, seq_error_rate, etc.)
///                  required to perform simulations.
/// - `pop`        : Optional tag defining the population of origin of the founder individual
#[derive(Debug, Clone)]
pub struct Pedigree {
    pub individuals: BTreeMap<String, Rc<RefCell<Individual>>>,
    pub comparisons: PedComparisons,
    params: Option<PedigreeParams>,
    pop : Option<String>,
}

impl Pedigree {
    /// Instantiate a blank pedigree.
    /// # @TODO -> Convert this to default();
    pub fn new() -> Pedigree {
        Pedigree { individuals: BTreeMap::new(), comparisons: PedComparisons::new(), params: None, pop: None}
    }

    /// Compute the pairwise differences of all the contained `self.comparisons`.
    /// # Arguments
    /// - `cont_af`            : Contaminating population allele frequencies for each compaired individual, at the current position.
    /// - `pileup_error_probs` : Sequencing error probabilities for the current position. Computed using the pileup phred-scale.
    /// 
    /// # Errors:
    /// - if `self.params` is `None`
    /// - if an individual being compared carries `None` alleles.
    #[inline]
    pub fn compare_alleles(&mut self, cont_af: [f64; 2], pileup_error_probs: &[f64; 2], rng: &mut fastrand::Rng) -> Result<()> {
        use PedigreeError::FailedAlleleComparison;
        // ---- Extract the user-defined contamination rate of this pedigree.
        let param = || self.get_params().with_loc(|| FailedAlleleComparison);
        let contam_rate = param()?.contam_rate;

        // ---- A 'None' pedigree param error rate implies the user did not provide any custom error_rate and/or 
        //      wishes to use the pileup local error rate.
        let seq_error_rate = match param()?.seq_error_rate {
            Some(error_rate) => error_rate,
            None             => *pileup_error_probs,
        };

        // ---- update the PWD of all comparisons at the current position.
        for comparison in &mut self.comparisons.iter_mut() {
            comparison.compare_alleles(contam_rate, cont_af, seq_error_rate, rng).with_loc(|| FailedAlleleComparison)?;
        }
        Ok(())
    }

    /// Populate the genotypes of all the founder alleles at the current positions.
    /// # Arguments:
    /// - `reader`: a `GenotypeReader`, either `VCFReader` or `FSTReader`
    /// - `rng`: random number generator. used to decide whether or not allele_frequency downsampling should be performed.
    /// 
    /// # Errors
    /// - if `self.params` is `None`
    /// 
    /// # Panics
    /// - if any founder's tag is set to `None`
    #[inline]
    pub fn update_founder_alleles(&mut self, reader: &dyn GenotypeReader, rng: &mut fastrand::Rng) -> Result<()> {
        let loc_msg = "While updating founder individiuals' alleles";
        // ---- Extract this pedigree allele frequency downsampling rate.
        let af_downsampling_rate = self.get_params().loc(loc_msg)?.af_downsampling_rate;

        // ---- Perform allele fixation at random, according to this pedigrees af_downsampling_rate
        if rng.f64() < af_downsampling_rate {
            self.founders_mut().for_each(|mut founder| founder.alleles = Some([0, 0]))
        } else {
            // ---- Fetch and assign the 'true' alleles for all founder individuals. 
            for mut founder in self.founders_mut() {

                // ---- Extract founder tag ; raise an error if None is returned.
                let founder_tag = founder.get_tag().loc(loc_msg)?;

                // ---- Fetch and assign the individual's allele from our reader.
                founder.alleles = Some(reader.get_alleles(founder_tag)?)
            }
        }
        Ok(())
    }

    /// Populate the genotypes of all the offspring alleles at the current position.
    /// # Arguments:
    /// - `inverval_prob_recomb`: probability that a recombination occured between the current and the previously typed position
    ///                           Computed through the use of a `GeneticMap` interval tree.
    /// - `pedigree_index`      : index of the current pedigree within the `PedigreeRep` vector. 
    ///                           Serves no purpose appart from logging and debugging
    /// 
    /// # Errors
    /// - if any individual's `self.strands` field is set tot `None`.
    /// 
    /// # Panics
    /// - if any individual's `self.parents` is set to none `None` (i.e. Individual is a founder.)
    #[inline]
    pub fn compute_offspring_alleles(&mut self, interval_prob_recomb: f64, pedigree_index: usize, rng: &fastrand::Rng, xchr_mode: bool) -> Result<()> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_alleles(interval_prob_recomb, pedigree_index, rng, xchr_mode)
            .with_loc(|| format!("While attempting to assign the alleles of {}", offspring.label))?;
        }
        Ok(())
    }

    ///  Wrap multiple simulations parameters within a new `PedigreeParam` struct and update `self.params` with it.
    pub fn set_params(&mut self, snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: Option<[f64; 2]>, contam_rate: [f64; 2]) {
        //trace!("error_rate: {seq_error_rate} | contam_rate: {contam_rate}");
        self.params = Some(
            PedigreeParams::new(snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate)
        );
    }

    /// Access this pedigree's parameters set.
    /// # Errors:
    /// - if `self.params` is `None`
    pub fn get_params(&self) -> Result<&PedigreeParams> {
        self.params.as_ref().with_loc(|| PedigreeError::EmptyParam)
    }

    /// Randomly assign the chromosome provenance of each offspring.
    /// # Errors:
    /// - if any of the offspring `self.parents` field is set to `None`
    pub fn assign_offspring_strands(&mut self) -> Result<()> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_strands().with_loc(||PedigreeError::FailedAlleleAssignment(offspring.label.clone()))?;
        }
        Ok(())
    }

    /// Randomly assign the sex of each individual.
    pub fn assign_random_sex(&mut self) -> Result<()> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_random_sex().with_loc(||PedigreeError::FailedSexAssignment(offspring.label.clone()))?;
        }
        Ok(())
    }

    /// Instantiate and include a new individual within this pedigree. 
    /// # Fields:
    /// - `label`  : name of the individual (e.g. "child")
    /// - `parents`: Optional name of the parents (e.g. ("mother", "father")). 
    ///              This should be set to `None` if the individual is a founder.
    ///
    /// # Errors: 
    /// - returns `std::io::result::InvalidInput` If any of the parents cannot be found within `self.individuals`
    pub fn add_individual(&mut self, label: &str, parents: Option<(&str, &str)>, sex: Option<Sex>) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parents = match parents {
            None => None,
            Some((parent1, parent2)) => {
                let parent1 = self.individuals.get(parent1).ok_or(InvalidInput)?;
                let parent2 = self.individuals.get(parent2).ok_or(InvalidInput)?;
                Some([parent1, parent2])
            },
        };
        let ind = Rc::new(RefCell::new(Individual::new(label, parents, sex)));
        self.individuals.insert(label.to_owned(), ind);
        Ok(())
    }

    /// Instantiate and include a new PedComparison within this pedigree.
    /// 
    /// # Fields
    /// - `label`: name of the comparison (e.g.: "Siblings", "Cousins", etc.)
    /// - `pair` : name of the individuals being compared (e.g.: ('daughter', 'son') ).
    /// 
    /// # Errors
    /// - returns `std::io::result::InvalidInput` If any individual used for the comparison can't 
    ///   be found within `self.individuals`
    pub fn add_comparison(&mut self, label: &str, pair: (&str, &str)) -> std::io::Result<()> {
        use std::io::ErrorKind::InvalidInput;
        let pair0 = self.individuals.get(pair.0).ok_or(InvalidInput)?;
        let pair1 = self.individuals.get(pair.1).ok_or(InvalidInput)?;
        self.comparisons.push(PedComparison::new(label, [pair0, pair1], pair0==pair1));
        Ok(())
    }

    /// Define the parents of a given offspring individual.
    /// # Arguments:
    /// - `ind`: name of the target individual.
    /// - `parents`: names of the target individual's parents.
    /// 
    /// # Errors 
    /// - returns `std::io::result::InvalidInput` If the target individual, or any of the parents cannot be found 
    ///   within `self.individuals`
    pub fn set_relationship(&mut self, ind: &str, parents: (&str, &str)) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parent0 = &self.individuals.get(parents.0).ok_or(InvalidInput)?.clone();
        let parent1 = &self.individuals.get(parents.1).ok_or(InvalidInput)?.clone();

        self.individuals.get_mut(ind)
            .ok_or(InvalidInput)?
            .borrow_mut()
            .set_parents([parent0, parent1]);
            Ok(())

    }

    /// Return a mutable borrow to a target individual. Only used during tests to bypass methods and mock individuals.
    /// # Arguments:
    /// - `label`: name of the target individual;
    #[cfg(test)]
    pub fn get_mutind(&mut self, label: &String) -> Result<std::cell::RefMut<Individual>, std::io::Error>{
        use std::io::ErrorKind::InvalidInput;
        Ok(self.individuals.get_mut(label)
            .ok_or(InvalidInput)?
            .borrow_mut())
    }

    /// Obtain a vector of mutable references leading to the founder individuals of this pedigree.
    /// 
    /// # @TODO:
    ///  - founders_mut and offsprings_mut have the same code, except the predicate:
    ///    - founters_mut    : ind.is_founder()
    ///    - offsprings_mut : !ind.is_founder()
    ///  - collecting the iterator is useless, since we're always calling this method to iterate over the items.
    ///    --> return an iterator!!
    #[inline]
    pub fn founders_mut(&mut self) -> impl Iterator<Item = RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
    }

    /// Obtain a vector of mutable references leading to the offsprings of this pedigree.
    #[inline]
    pub fn offsprings_mut(&mut self) -> impl Iterator<Item = RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| !RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
    }

    /// Set the population tags, and assign random SampleTag for each founder individual within this pedigree.
    /// # Arguments:
    /// - `panel`: input samples panel definition, in the form of a `PanelReader`.
    /// - `pop`  : (super-)population id requested for this pedigree's founder individuals.
    /// - `contaminants`: Set of contaminating individuals, in the form of a `Contaminant`
    /// 
    /// #Errors:
    /// - if `pop` name is invalid.
    pub fn set_founder_tags(&mut self, panel: &PanelReader, pop: &String, contaminants: Option<&Contaminant>) -> Result<()>{
        use PedigreeError::{MissingContaminant, MissingSampleTag};
        self.pop = Some(pop.to_owned());

        // ---- Contaminating individual are excluded from pedigree individuals.
        let mut exclude_tags = contaminants.map(|cont| cont.as_flat_list()).with_loc(||MissingContaminant)?;

        // ---- For each founder, pick and assign a random SampleTag using our panel (without replacement)
        for mut founder in self.founders_mut() {

            // ---- Pick a random SampleTag from our panel.
            let random_tag = panel.random_sample(pop, Some(&exclude_tags), founder.sex)?.with_loc(||MissingSampleTag)?; 
            founder.set_tag(random_tag.clone());
            exclude_tags.push(random_tag); // Exclude this pedigree individual for other iterations.
        };
        Ok(())
    }

    /// Set all individual allele's to `None` within this pedigree.
    pub fn clear_alleles(&mut self){
        for ind in self.individuals.values_mut(){
            ind.borrow_mut().clear_alleles()
        }
    }

    pub fn all_sex_assigned(&self) -> bool {
        self.individuals.values().all(|ind| ind.borrow().is_sex_assigned())
    }

}

impl Default for Pedigree {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn test_pedigree_set() -> std::io::Result<Pedigree> {
        let mut pedigree = Pedigree::new();
        pedigree.add_individual("father", None, None)?;
        pedigree.add_individual("mother", None, None)?;
        pedigree.add_individual("offspr", Some(("father", "mother")), None)?;
    
        let mut father = pedigree.get_mutind(&"father".to_string()).expect("Cannot extract father");
        father.set_alleles([0, 1]);
        drop(father);
    
        let mut mother = pedigree.get_mutind(&"mother".to_string()).expect("Cannot extract mother");
        mother.set_alleles([1, 0]);
        drop(mother);
    
    
        Ok(pedigree)
    }

    type TestPedDef = Vec<(&'static str, Option<(&'static str, &'static str)>)>;
    
    fn test_pedigree_random(def: Option<TestPedDef>) -> std::io::Result<Pedigree> {
        let mut pedigree = Pedigree::new();
        if let Some(map) = def {
            for (label, parents) in map {
                println!("{label}, {parents:?}");
                pedigree.add_individual(label, parents, None)?;
            }
            for (_label, ind) in pedigree.individuals.iter_mut() {
                if ind.borrow().is_founder() {
                    ind.borrow_mut().set_alleles([fastrand::bool() as u8, fastrand::bool() as u8]);
                }
            }
        } else {
            let def = Vec::from([
                ("father", None),
                ("mother", None),
                ("offspr", Some(("father", "mother")))
            ]);
            return test_pedigree_random(Some(def))
        }

        Ok(pedigree)
    }


    #[test]
    #[should_panic]
    fn meiosis_assign_alleles_empty_strands(){
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.assign_alleles(0.0, 0, &fastrand::Rng::new(), false).expect("Failed to assign alleles");
    }

    #[test]
    fn meiosis_assign_alleles_filled_strands(){
        let rng          = fastrand::Rng::new();
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        let mut offspr   = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands   = Some([0, 0]);

        let output = offspr.assign_alleles(0.0, 0, &rng, false).expect("Failed to assign alleles");
        assert!(output);

        let output = offspr.assign_alleles(0.0, 0, &rng, false).expect("Failed to assign alleles");
        assert!(!output);
    }

    #[test]
    fn meiosis_check_strands_00() {
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        let mut offspr   = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands   = Some([0, 0]);
        offspr.assign_alleles(0.0, 0, &fastrand::Rng::new(), false).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([0, 1]))
    }

    #[test]
    fn meiosis_check_strands_01() {
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        let mut offspr   = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands   = Some([1, 1]);

        offspr.assign_alleles(0.0, 0, &fastrand::Rng::new(), false).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 0]));
        assert_eq!(offspr.currently_recombining, [false, false]);

    }

    #[test]
    fn meiosis_check_recombination() {
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        let mut offspr   = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands   = Some([0, 1]);
        offspr.assign_alleles(1.0, 0, &fastrand::Rng::new(), false).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 1]));
        assert_eq!(offspr.currently_recombining, [true, true]);

    }

    #[test]
    fn assign_offspring_strands() -> Result<()> {
        let mut pedigree = test_pedigree_set().expect("Cannot generate test pedigree");
        pedigree.assign_offspring_strands()?;
        Ok(())
    }

    #[test]
    fn assign_sex() -> Result<()> {
        let def = Vec::from([
            ("F1.1", None),
            ("F1.2", None),
            ("F1.3", None),
            ("F2.1", None),
            ("O2.2", Some(("F1.1", "F1.2"))),
            ("O2.3", Some(("F1.1", "F1.2"))),
            ("F2.4", None),
            ("02.5", Some(("F1.2", "F1.3"))),
            ("O3.1", Some(("F2.1", "O2.2"))),
            ("O3.2", Some(("O2.3", "F2.4")))
            
        ]);
        for _ in 0..1000 {
            let mut pedigree = test_pedigree_random(Some(def.clone())).unwrap();//.expect("Cannot generate test pedigree");
            pedigree.assign_random_sex()?;
            for (_label, ind) in pedigree.individuals.iter() {
                assert!(ind.borrow().sex.is_some());
            }
        }
        Ok(())
    }
}