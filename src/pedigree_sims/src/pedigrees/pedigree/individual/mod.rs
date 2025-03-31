use std::{
    cell::RefCell,
    cmp::{Ord, Ordering, PartialOrd},
    hash::{Hash, Hasher},
    rc::Rc,
};

use grups_io::read::SampleTag;
use genome::Sex;
use fastrand;

use located_error::prelude::*;
use log::trace;

mod parents;
use parents::Parents;

use self::error::IndividualError;

mod error;

/// Space padding lengths used for `std::fmt::Display` of Individual
const TAG_DISPLAY_LEN    : usize = 10; // Space padding of `self.tag`
const LABEL_DISPLAY_LEN  : usize = 10; // Space padding of `self.label`
const PARENTS_DISPLAY_LEN: usize = 25; // Space padding of `self.parents`

/// Pedigree Individual.
/// # Fields:
/// - `tag`                  : Optional SampleTag of the Individual 
///                              - `Some(SampleTag)` if the individual is a founder.
///                              - `None`            if the individual is a simulated offspring.
/// - `label`                : User-defined name of the individual (e.g. 'child', 'father', 'mother')
/// - `parents`              : Optional Set of references for each parents of the individual.
///                              - `None`             if the individual is a founder.
///                              - `Some(references)` if the individual is a simulated offspring.
/// - `strands`              : Optional array indicating the provenance of each individual's chromosome strand.
///                            Strand provenance is tracked as a set of two indices. `strands[i]` = strand index of 
///                            `parents[i]`. Two possible values: `0` = left strand and `1` == right strand of the 
///                            given parent. 
///                              - `None`          if the individual is a founder.
///                              - `Some(strands)` if the individual is a simulated offspring. 
/// - `currently_recombining`: Array of `bool` tracking whether or not the parent's chromosome is currently
///                            recombining. `currently_recombining[i]` = tracker for `parents[i]`.
/// - `alleles`              : Optional size-two set of alleles of the Individual for the current SNP position.
#[derive(Debug, Clone)]
pub struct Individual {
    pub tag                  : Option<SampleTag>,
    pub label                : String,
    parents                  : Option<Parents>,
    pub strands              : Option<[usize; 2]>,
    pub currently_recombining: [bool; 2],
    pub alleles              : Option<[u8; 2]>,
    pub sex                  : Option<Sex>,
}

impl std::fmt::Display for Individual {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let parents = match &self.parents {
            None => "None".to_string(),
            Some(parents) => format!("{parents}")
        };
        let tag = match &self.tag {
            Some(tag) => tag.id().clone(),
            None => "None".to_owned()
        };
        write!(f, "tag: {: <TAG_DISPLAY_LEN$} label: {: <LABEL_DISPLAY_LEN$} - parents: {: <PARENTS_DISPLAY_LEN$}", tag, self.label, parents)
    }
}

impl PartialEq for Individual {
    fn eq(&self, other: &Individual) -> bool {
        self.label == other.label
    }
}

impl Eq for Individual {}

impl Hash for Individual {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.label.hash(state);
    }
}

impl std::borrow::Borrow<String> for Individual {
    fn borrow(&self) -> &String {
        self.label.borrow()
    }
}

impl Ord for Individual {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.label).cmp(&(other.label))
    }
}

impl PartialOrd for Individual {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

type ParentsRef<'a> = [&'a Rc<RefCell<Individual>>; 2];

impl Individual {
    /// Instantiate a new individual.
    /// # Arguments
    /// - `label`  : User-defined name of the individual (e.g. "father", "mother", "child", etc.)
    /// - `parents`: Size-two array of `&Rc<RefCell<Individual>>`, representing the individual's parents. 
    pub fn new(label: &str, parents: Option<ParentsRef>, sex: Option<Sex> ) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label: label.to_string(), parents, strands: None, currently_recombining: [false, false], alleles: None, sex}
    }

    /// Manually set the alleles for this individual. Used during tests only, to bypass some methods and create mock Individuals.
    #[cfg(test)]
    pub fn set_alleles(&mut self, alleles: [u8; 2]) {
        self.alleles = Some(alleles);
    }

    /// Return the contents of `self.alleles`
    /// # Errors
    /// - In case `self.alleles` is `None`
    pub fn get_alleles(&self) -> Result<[u8; 2]> {
        self.alleles.with_loc(||IndividualError::MissingAlleles)
    }

    /// Simulate meiosis for the current position and return a unique allele (as `u8`).
    /// Arguments:
    /// - `selected_strand`                : index of `self.alleles` used for sampling. The provided value will most likely
    ///                                      originate from `self.strands` of the simulated offspring. 
    /// - `offspring_currently_recombining`: bool indicating whether the strand being currently simulated is under recombination.
    ///                                      The provided value will most likely originate from `self.currently_recombining` of
    ///                                      the simulated offspring.
    /// # Panics:
    /// - when `self.alleles` is `None`
    pub fn meiosis(&self, selected_strand: usize, offspring_currently_recombining: bool) -> u8 {
        // --- Switch the selected strand if the offspring is currently recombining.
        let selected_strand = match offspring_currently_recombining {
            false => selected_strand,
            true => (selected_strand + 1) % 2
        };

        // ---- Return the allele at index = selected_strand. 
        match self.alleles {
            None          => panic!("Trying to perform meiosis within an empty genome!"),
            Some(alleles) => alleles[selected_strand],
        }
    }

    /// Randomly assign a strand provenance for each parents. This method should be called once, before the simulations.
    /// Returns `true` if strand assignment was succesful and non-redundant, `false` if the individual's strand were 
    /// already assigned.
    /// 
    /// # Errors:
    /// - if `self.parents` is `None`
    pub fn assign_strands(&mut self) -> Result<bool> {
        if self.parents.is_none() {
            return loc!(IndividualError::MissingParents)
        }

        if self.strands.is_some() { // Strands are already assigned. Skip.
            return Ok(false)
        }
        let rng = fastrand::Rng::new();
        self.strands = Some([rng.usize(0..=1), rng.usize(0..=1)]);
        Ok(true)
    }

    /// Return a reference to `self.tag`
    pub fn get_tag(&self) -> Option<&SampleTag> {
        self.tag.as_ref()
    }

    /// Manually set the Individuals parents.
    /// # Arguments
    /// - `parents`: Size-two array of `&Rc<RefCell<Individual>>`, representing the individual's parents. 
    pub fn set_parents(&mut self, parents: ParentsRef) {
        self.parents = Some(Self::format_parents(parents));
    }

    /// Rc::clone() the provided pair of of parents during instantiation. (see. `Individual::new()`)
    fn format_parents(parents:  [&Rc<RefCell<Individual>>; 2]) -> Parents {
        Parents::new([Rc::clone(parents[0]), Rc::clone(parents[1])])
    }

    /// Check whether or not this individual is a founder individual. Returns `true` if `self.parents == None`
    #[inline(always)]
    pub fn is_founder(&self) -> bool {
        self.parents.is_none()
    }

    /// Manually set the individuals SampleTag. meaningless for offsprings.*
    /// Arguments:
    /// - `tag`: Input sample identification tag. (i.e. the name and idx of the snp-callset individual being used for simulations)
    pub fn set_tag(&mut self, tag: SampleTag){
        self.tag = Some(tag);
    }

    /// Set the Individual's `self.alleles` field to `None`. 
    /// This method is most likely called when switching from one simulated SNP coordinate to another.
    pub fn clear_alleles(&mut self){
        self.alleles = None
    }
    
    /// Set the individual's alleles, given a recombination probability. 
    /// - Returns `true`  if assignment was successful and non-redundant.
    /// - Returns `false` if `self.alleles` were already assigned.
    /// 
    /// # Arguments:
    /// - `recombination_prob`: probability that a recombination occured between the previous and current coordinate.
    ///                         This was most likely computed from a genetic_map. See `genome::GeneticMap` and 
    ///                         `genome::RecombinationRange`
    /// - `ped_idx`           : replicate index of this individual's belonging pedigree. This parameter is purely provided
    ///                         for debugging and logging, and serves no purpose during allele assignment.
    /// # Errors
    /// - when trying to assign alleles while `self.strands` is `None`.
    /// 
    /// # Panics
    /// - When attempting to assign alleles while `self.parents` is `None` (i.e. Individual is a founder.)
    /// 
    /// # @ TODO
    /// - Instantiating a new Rng for each individual might not be very efficient...
    ///   Passing a &ThreadRng reference around might be better.
    pub fn assign_alleles(&mut self, recombination_prob: f64, ped_idx: usize, rng: &fastrand::Rng, xchr_mode: bool) -> Result<bool> {
        use IndividualError::{InvalidAlleleAssignment, MissingParents, MissingStrands};
        // ---- Ensure this method call is non-redundant.
        if self.alleles.is_some() {
            return Ok(false)
        }
        // ---- Ensure the individual is 'equipped' with parents and strands before attempting allele assignment.
        let Some(parents) = &self.parents else {
            return Err(anyhow!(MissingParents)).with_loc(||InvalidAlleleAssignment)
        };

        // ---- Perform allele assignment for each parent.
        for (i, parent) in parents.iter().enumerate() {
            // ---- Assign parent genome if not previously generated.
            if parent.borrow().alleles.is_none() {
                parent.borrow_mut().assign_alleles(recombination_prob, i, rng, xchr_mode)
                    .with_loc(||InvalidAlleleAssignment)?;
            }

            // ---- Check if recombination occured for each parent and update recombination tracker if so.
            if (!xchr_mode || parent.borrow().sex == Some(Sex::Female)) && rng.f64() < recombination_prob { 
                trace!("- Cross-over occured in ped: {:<5} - ind: {} ({} {:?})", ped_idx, self.label, parent.borrow().label, parent.borrow().sex);
                self.currently_recombining[i] = ! self.currently_recombining[i];
            }
        }
        
        // ---- Perform allele assignment for `self`, by simulating meiosis for each parent.
        let Some(strands) = self.strands else {
            return Err(anyhow!(MissingStrands)).with_loc(||InvalidAlleleAssignment)
        };

        self.alleles = if xchr_mode {
            let mut alleles = [0u8 ; 2];
            // ---- Find the index of both parents
            let father_idx = parents.iter().position(|p| p.borrow().sex == Some(Sex::Male)).expect("No parent found..");
            let mother_idx = (father_idx + 1 ) % 2;

            alleles[mother_idx] = parents[mother_idx].borrow().meiosis(strands[mother_idx], self.currently_recombining[mother_idx]);
            alleles[father_idx] = match self.sex {
                Some(Sex::Male)           => Ok(alleles[mother_idx]), // If the descendant is a male, alleles are exclusively from the mother
                Some(Sex::Female)         => Ok(parents[father_idx].borrow().meiosis(strands[father_idx], self.currently_recombining[father_idx])),
                Some(Sex::Unknown) | None => Err(IndividualError::UnknownOrMissingSex).loc("While attempting to assign alleles during X-chromosome-mode"),
            }?;

            // ---- Sanity checks.
            if self.sex == Some(Sex::Male) && alleles[0] != alleles[1] {
                // Male individuals are not expected to be heterozygous during X-chromosome simulations.
                return Err(IndividualError::SpuriousAlleleAssignment{alleles})
                    .loc("Male Individual is heterozygous while in X-chromosome mode")
            }
            if self.currently_recombining[father_idx] {
                // Fathers are not expected to recombine during X-chromosome simulations.
                return Err(IndividualError::InvalidOrSpuriousRecombinationEvent)
                    .loc("Father is recombining while in X-chromosome-mode")
            }

            Some(alleles)
        } else {
            let haplo_0 = parents[0].borrow().meiosis(strands[0], self.currently_recombining[0]);
            let haplo_1 = parents[1].borrow().meiosis(strands[1], self.currently_recombining[1]);
            Some([haplo_0, haplo_1])
        };

        Ok(true)
    }

    pub fn assign_random_sex(&mut self) -> Result<bool> {
        use IndividualError::InvalidSexAssignment;
        // ---- Ensure this method call is non-redundant
        if self.sex.is_some() {
            return Ok(false)
        }

        // ---- Perform sex-asignment for each parent (if the individual has known parents.)
        if let Some(parents) = &self.parents { 
            for (i, parent) in parents.iter().enumerate() {
                // ---- Assign parent sex if not previously decided.
                let spouse = &parents[(i+1) % 2];
                if parent.borrow().sex.is_none() { // If the parent's sex is still unknown
                    parent.borrow_mut().sex = if let Some(spouse_sex) = spouse.borrow().sex { // If the spouse sex is already known, assign the opposite sex.
                        match spouse_sex {
                            Sex::Female  => Some(Sex::Male),
                            Sex::Male    => Some(Sex::Female),
                            Sex::Unknown => None
                        }
                    } else {
                        Some(Sex::random()) // If not, assign a random sex to the parent.
                    }
                }
                // ---- Apply the same process for the parent.
                parent.borrow_mut().assign_random_sex()
                    .with_loc(||InvalidSexAssignment)?;
            }
        }

        // ---- Randomly assign sex of the considered individual
        self.sex = Some(Sex::random());
        Ok(true)
    }

    pub fn is_sex_assigned(&self) -> bool {
        self.sex.is_some_and(|s| s != Sex::Unknown)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pedigrees::pedigree::tests::common;

    fn perform_allele_asignment(offspring: &mut Individual, parents_alleles: [[u8;2];2], recombination_prob: f64) -> Result<()> {
        let rng = fastrand::Rng::new();
        let parents = offspring.parents.as_ref().expect("Missing parents");
        parents[0].borrow_mut().alleles = Some(parents_alleles[0]);
        parents[1].borrow_mut().alleles = Some(parents_alleles[1]);        
        offspring.assign_alleles(recombination_prob, 0, &rng, false)?;
        Ok(())
    }

    fn run_all_allele_assignment_cases(recombination_prob: f64) -> Result<()> {
        let mut offspring = common::mock_offspring("offspring", None);
        let valid_alleles = [[0,0], [0,1], [1,0], [1,1]];
        let mut valid_strands = [[0,0], [0,1], [1,0], [1,1]];

        for parent_0_alleles in valid_alleles.iter(){
            for parent_1_alleles in valid_alleles.iter() {
                for strands in valid_strands.iter_mut() {
                    offspring.strands = Some(*strands);
                    perform_allele_asignment(&mut offspring, [*parent_0_alleles, *parent_1_alleles], recombination_prob)?;

                    // If the individual's parent is 'recombining', we expect strand assignment to be inverted. 0 becomes 1 ; 1 becomes 0
                    for i in [0,1]{
                        if offspring.currently_recombining[i] {
                            strands[i] = (strands[i]+1) % 2;
                        }
                    }

                    let got = offspring.alleles.expect("Missing alleles within offspring.");
                    let want = [parent_0_alleles[strands[0]], parent_1_alleles[strands[1]]];

                    assert_eq!(got, want);
                    offspring.alleles = None;
                }
            }
        }
        Ok(())
    }

     #[test]
    fn alleles_getter_filled() -> Result<()> {
        let mut ind = common::mock_founder("offspring");
        let alt_ref = [0,1];
        for i in alt_ref {
            for j in alt_ref {
               let alleles = [i, j];
               ind.set_alleles(alleles);
               assert_eq!(ind.get_alleles()?, alleles);
            }
        }
        Ok(())
    }

    #[test]
    fn alleles_getter_empty(){
        let ind = common::mock_founder("offspring");
        let alleles = ind.get_alleles();
        assert!(alleles.is_err());
    }

    #[test]
    fn meiosis_not_recombining() {
        let offspring_currently_recombining = false; 
        let mut ind = common::mock_founder("offspring");
        let alt_ref = [0,1];
        for i in alt_ref {
            for j in alt_ref {
                for selected_strand in [0, 1] {

                    let alleles = [i, j];
                    ind.set_alleles(alleles);

                    let want = alleles[selected_strand];
                    let got = ind.meiosis(selected_strand, offspring_currently_recombining);
                    assert_eq!(want, got);
                }
            }
        }
    }

    #[test]
    fn meiosis_recombining() {
        let offspring_currently_recombining = true; 
        let mut ind = common::mock_founder("parent");
        let alt_ref = [0,1];
        for i in alt_ref {
            for j in alt_ref {
                for selected_strand in [0, 1] {
                    
                    let mut alleles = [i, j];
                    ind.set_alleles(alleles);

                    alleles.reverse();
                    let want = alleles[selected_strand];
                    let got = ind.meiosis(selected_strand, offspring_currently_recombining);
                    assert_eq!(want, got);
                }
            }
        }
    }

    #[test]
    #[should_panic]
    fn meiosis_empty_alleles() {
        let ind = common::mock_founder("parent");
        ind.meiosis(0, false);
    }


    #[test]
    fn strand_setter_offspring() -> Result<()> {
        let valid_strands = [[0,0], [0,1], [1,0], [1,1]];
        let mut ind = common::mock_offspring("offspring", None);

        for _ in 0..1000 {
            ind.assign_strands()?;
            let rng_strands = ind.strands.expect("Missing individual strands.");
            assert!(valid_strands.contains(&rng_strands));
        }
        Ok(())
    }

    #[test]
    fn strand_setter_founder() {
        let mut ind = common::mock_founder("parent");
        let result = ind.assign_strands();
        assert!(result.is_err())
    }

    #[test]
    fn sex_setter_founder() -> Result<()>{
        let mut ind = common::mock_founder("parent");
        ind.assign_random_sex()?;
        assert!(ind.sex.is_some());
        Ok(())
    }

    #[test]
    fn sex_setter_offspring() -> Result<()> {
        for _ in 0..1000 {
            let mut child = common::mock_offspring("child", Some(["parent-1", "parent-2"]));
            child.assign_random_sex()?;
            let parents = child.parents.expect("Test offspring has missing parents");
            assert_ne!(parents[0].borrow().sex, parents[1].borrow().sex);
        }
        Ok(())
    }

    #[test]
    fn get_set_sampletag() {
        let mut ind = common::mock_founder("parent");
        let tag = SampleTag::new("HG00096", Some(0), None);
        ind.set_tag(tag.to_owned());
        assert_eq!(ind.get_tag(), Some(&tag));
        
    }

    #[test]
    fn get_empty_tag(){
       let ind = common::mock_founder("parent");
       let result = ind.get_tag();
       assert!(result.is_none());
    }

    #[test]
    fn founder_is_founder() {
        let ind = common::mock_founder("parent");
        assert!(ind.is_founder());
    }

    #[test]
    fn offpsring_is_not_founder() {
        let ind = common::mock_offspring("offspring", None);
        assert!(!ind.is_founder());
    }

    #[test]
    fn clear_alleles() {
        let mut ind = common::mock_founder("parent");
        ind.set_alleles([0,1]);
        assert!(ind.alleles.is_some());
        ind.clear_alleles();
        assert!(ind.alleles.is_none());


    }

    #[test]

    fn alleles_assignment_founder() {
        let mut ind = common::mock_founder("parent");
        let result = ind.assign_alleles(0.0, 0, &fastrand::Rng::new(), false);
        assert!(result.is_err())
    }

    #[test]
    fn alleles_assignments_unnassigned_parent_alleles(){
        let mut ind = common::mock_offspring("offspring", None);
        let result = ind.assign_alleles(0.0, 0, &fastrand::Rng::new(), false);
        assert!(result.is_err())
    }

    #[test]
    fn alleles_assignment_offspring_no_recombination() -> Result<()> {
        run_all_allele_assignment_cases(0.0)
    }
    #[test]
    fn alleles_assignment_offspring_full_recombination() -> Result<()> {
        run_all_allele_assignment_cases(1.0)
    }

    #[test]
    fn allele_assignment_updates_recombination_status() -> Result<()> {
        let mut offspring = common::mock_offspring("offspring", None);
        offspring.strands = Some([0,0]);
        let parents_alleles = [[0,1], [0,1]];
        let recombination_prob = 1.0;

        assert_eq!(offspring.currently_recombining, [false, false]);
        perform_allele_asignment(&mut offspring, parents_alleles, recombination_prob)?;
        assert_eq!(offspring.currently_recombining, [true, true]);
        Ok(())
    }

    #[test]
    fn ind_equality() {
        let  ind1 = common::mock_founder("parent");
        let  ind2 = common::mock_founder("parent");
        assert_eq!(ind1, ind2)
    }

    #[test]
    fn ind_inequality() {
        let  ind1 = common::mock_founder("ind1");
        let  ind2 = common::mock_founder("ind2");
        assert_ne!(ind1, ind2)
    }

    #[test]
    fn hashable() {
        let mut ind_set = std::collections::HashSet::new();
        let n_iters: u32 = 10_000;
        for i in 0..n_iters {
            let  new_ind = common::mock_founder(&i.to_string());
            assert!(ind_set.insert(new_ind.clone()));
            assert!(ind_set.contains(&new_ind));
        }
    }

    #[test]
    fn ordering() {
        let  ind_a = common::mock_founder("A");
        let  ind_b = common::mock_founder("B");
        assert!(ind_a <  ind_b);
        assert!(ind_b <= ind_b);

        assert!(ind_b >= ind_a);
        assert!(ind_b >  ind_a);

        // Different Index should not impact ordering. What matters is the ID.
        let mut ind_a_prime = ind_a.clone();
        ind_a_prime.set_tag(SampleTag::new("A", Some(996), None));
        assert!(ind_a <= ind_a_prime);


    }

    #[test]
    fn display() {
        let (offspring_label, father_label, mother_label) = ("ind1", "ind2", "ind3");
        let mut offspring = common::mock_offspring(offspring_label, Some([father_label, mother_label]));
        
        let display = format!("{offspring}");
        assert!(display.contains(offspring_label));
        assert!(display.contains(father_label));
        assert!(display.contains(mother_label));

        offspring.set_tag(SampleTag::new("HG00096", None, None));
        let display = format!("{offspring}");
        assert!(display.contains("HG00096"));
    }
}