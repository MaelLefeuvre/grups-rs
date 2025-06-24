use std::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt::{self, Display, Formatter},
};

use std::io::Write;


#[cfg(test)] 
use std::hash::{Hash, Hasher};

use grups_io::read::SampleTag;
use genome::Sex;
use fastrand;

use located_error::prelude::*;
use log::trace;

mod error;
pub use error::IndividualError;
use slotmap::KeyData;

use super::RelationshipId;

//mod error;
/// Space padding lengths used for `std::fmt::Display` of Individual
const TAG_DISPLAY_LEN    : usize = 10; // Space padding of `self.tag`
const LABEL_DISPLAY_LEN  : usize = 10; // Space padding of `self.label`
const PARENTS_DISPLAY_LEN: usize = 25; // Space padding of `self.parents`

// --------------------------------------------------------------------------------------------- //
// ---- Individual Id
slotmap::new_key_type! {pub struct IndividualId;}

impl IndividualId {
    pub fn to_u64(self) -> u64 {
        self.0.as_ffi()
    }

    pub fn from_u64(id: u64) -> Self {
        IndividualId::from(KeyData::from_ffi(id))
    }
}

impl Display for IndividualId {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.0)
    }
}

#[derive(Debug, Clone)]
pub struct Individual {
    pub id: IndividualId,
    tag: Option<SampleTag>,
    label: String,
    parents: Option<[RelationshipId; 2]>,
    pub strands: Option<[usize; 2]>,
    pub currently_recombining: [bool;2],
    pub alleles: Option<[u8; 2]>,
    pub sex: Option<Sex>
}

impl Display for Individual {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)
    }
}

impl Individual {
    pub fn new(id: IndividualId, label: &str, sex: Option<Sex>) -> Self {
        Self {
            id,
            tag: None,
            label: label.to_string(),
            parents: None,
            strands: None,
            currently_recombining: [false, false],
            alleles: None,
            sex
        }
    }

    #[cfg(test)]
    pub fn set_alleles(&mut self, alleles: [u8; 2]) {
        self.alleles = Some(alleles);
    }

    pub fn get_alleles(&self) -> Result<[u8; 2]> {
        self.alleles.with_loc(|| IndividualError::MissingAlleles)
    }

    pub fn add_parents(&mut self, relationships: [RelationshipId; 2]) {
        self.parents = Some(relationships);
    }
    pub fn has_parents(&self) -> bool {
        self.parents.is_some()
    }

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

    pub fn assign_strands(&mut self) -> Result<bool> {
        if self.parents.is_none() {
            return loc!(IndividualError::MissingParents)
        }

        if self.strands.is_some() { // Strands are already assigned. Skip.
            return Ok(false)
        }
        let mut rng = fastrand::Rng::new();
        self.strands = Some([rng.usize(0..=1), rng.usize(0..=1)]);
        Ok(true)
    }

    /// Return a reference to `self.tag`
    pub fn get_tag(&self) -> Option<&SampleTag> {
        self.tag.as_ref()
    }

    /// Return the indivdiual's label
    pub fn label(&self) -> &str {
        self.label.as_str()
    }

    /// Manually set the Individuals parents.
    /// # Arguments
    /// - `parents`: Size-two array of `&Arc<RwLock<Individual>>`, representing the individual's parents. 
    pub fn set_parents(&mut self, parents: [RelationshipId; 2]) {
        self.parents = Some(parents);
    }

    pub fn get_parents(&self) -> Option<[RelationshipId; 2]> {
        self.parents
    }

    /// Check whether or not this individual is a founder individual. Returns `true` if `self.parents == None`
    #[allow(clippy::inline_always)]
    #[inline(always)]
    pub fn is_founder(&self) -> bool {
        self.parents.is_none()
    }

    /// Check whether or not this individual is a founder individual. Returns `true` if `self.parents == None`
    #[allow(clippy::inline_always)]
    #[inline(always)]
    pub fn is_offspring(&self) -> bool {
        self.parents.is_some()
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
        self.alleles = None;
    }

    pub fn is_sex_assigned(&self) -> bool {
        self.sex.is_some_and(|s| s != Sex::Unknown)
    }

}



#[cfg(test)]
mod tests {
    use fastrand::Rng;

    use super::*;
    use crate::pedigrees::pedigree::tests::common;
    use crate::pedigrees::Pedigree;

    fn perform_allele_assignment(pedigree: &mut Pedigree, iid: IndividualId, parents_alleles: [[u8;2];2], recombination_prob: f64) -> Result<bool> {
        let mut rng = fastrand::Rng::new();
        
        let ind = pedigree.get_ind(iid).unwrap();
        let parent_rels = ind.get_parents().expect("Missing parents");
        for i in [0, 1] {
            let parent_id = pedigree.edges.get(parent_rels[i]).unwrap().to;
            pedigree.get_ind_mut(parent_id).unwrap().alleles =  Some(parents_alleles[i]);
        }
        pedigree.assign_alleles(iid, recombination_prob, 0, false, &mut rng)
    }

    fn run_all_allele_assignment_cases(recombination_prob: f64) -> Result<()> {
        let mut offspring_pedigree = common::mock_offspring_pedigree("offspring", None);
        let offspring_id = offspring_pedigree.get_ind_id("offspring").unwrap();

        let valid_alleles = [[0,0], [0,1], [1,0], [1,1]];
        let mut valid_strands = [[0,0], [0,1], [1,0], [1,1]];

        for parent_0_alleles in &valid_alleles{
            for parent_1_alleles in &valid_alleles {
                for strands in &mut valid_strands {

                    offspring_pedigree.get_ind_mut(offspring_id).unwrap().strands = Some(*strands);
                    perform_allele_assignment(&mut offspring_pedigree, offspring_id, [*parent_0_alleles, *parent_1_alleles], recombination_prob)?;

                    // If the individual's parent is 'recombining', we expect strand assignment to be inverted. 0 becomes 1 ; 1 becomes 0
                    let offspring = offspring_pedigree.get_ind(offspring_id).unwrap();
                    for (i, strand) in strands.iter_mut().enumerate() {
                        if offspring.currently_recombining[i] {
                            *strand = (*strand + 1) % 2;
                        }
                    }

                    let got  = offspring.alleles.expect("Missing alleles within offspring.");
                    let want = [parent_0_alleles[strands[0]], parent_1_alleles[strands[1]]];
                    println!("{got:?} | {want:?}");
                    assert_eq!(got, want);
                    offspring_pedigree.get_ind_mut(offspring_id).unwrap().alleles = None;
                }
            }
        }
        Ok(())
    }

     #[test]
    fn alleles_getter_filled() -> Result<()> {
        let mut pedigree = common::mock_founder_pedigree("offspring");
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
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
        let mut pedigree = common::mock_founder_pedigree("offspring");
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
        let alleles = ind.get_alleles();
        assert!(alleles.is_err());
    }

    #[test]
    fn meiosis_not_recombining() {
        let mut pedigree = common::mock_founder_pedigree("offspring");
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
        let offspring_currently_recombining = false; 
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
        let mut pedigree = common::mock_founder_pedigree("offspring");
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
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
    #[should_panic = "Trying to perform meiosis within an empty genome!"]
    fn meiosis_empty_alleles() {
        let mut pedigree = common::mock_founder_pedigree("offspring");
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
        ind.meiosis(0, false);
    }


    #[test]
    fn strand_setter_offspring() -> Result<()> {
        let valid_strands = [[0,0], [0,1], [1,0], [1,1]];
        let mut pedigree = common::mock_offspring_pedigree("offspring", None);
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();

        for _ in 0..1000 {
            ind.assign_strands()?;
            let rng_strands = ind.strands.expect("Missing individual strands.");
            assert!(valid_strands.contains(&rng_strands));
        }
        Ok(())
    }

    #[test]
    fn strand_setter_founder() {
        let mut pedigree = common::mock_founder_pedigree("parent");
        let ind = pedigree.get_ind_from_label_mut("parent").unwrap();
        let result = ind.assign_strands();
        assert!(result.is_err());
    }

    #[test]
    fn sex_setter_founder() -> Result<()>{
        let mut pedigree = common::mock_founder_pedigree("parent");
        let iid = pedigree.get_ind_from_label("parent").unwrap().id;
        pedigree.assign_random_sex(iid)?;
        assert!(pedigree.get_ind(iid).unwrap().sex.is_some());
        Ok(())
    }

    #[test]
    fn sex_setter_offspring() -> Result<()> {
        for _ in 0..1000 {
            let mut pedigree = common::mock_offspring_pedigree("child", Some(["parent-1", "parent-2"]));
            let child_id = pedigree.get_ind_from_label_mut("child").unwrap().id;
            pedigree.assign_random_sex(child_id)?;
            let parents = pedigree.get_parents(child_id).expect("Test offspring has missing parents");
    
            println!("hey: {parents:#?}");
            assert_ne!(parents[0].sex, parents[1].sex);
        }
        Ok(())
    }

    #[test]
    fn get_set_sampletag() {
        let mut pedigree = common::mock_founder_pedigree("parent");
        let ind = pedigree.get_ind_from_label_mut("parent").unwrap();
        let tag = SampleTag::new("HG00096", Some(0), None);
        ind.set_tag(tag.clone());
        assert_eq!(ind.get_tag(), Some(&tag));
        
    }

    #[test]
    fn get_empty_tag(){
        let mut pedigree = common::mock_founder_pedigree("parent");
        let ind = pedigree.get_ind_from_label_mut("parent").unwrap();
        let result = ind.get_tag();
        assert!(result.is_none());
    }

    #[test]
    fn founder_is_founder() {
        let mut pedigree = common::mock_founder_pedigree("parent");
        let ind = pedigree.get_ind_from_label_mut("parent").unwrap();
        assert!(ind.is_founder());
    }

    #[test]
    fn offspring_is_not_founder() {
        let mut pedigree = common::mock_offspring_pedigree("offspring", None);
        let ind = pedigree.get_ind_from_label_mut("offspring").unwrap();
        assert!(!ind.is_founder());
    }

    #[test]
    fn clear_alleles() {
        let mut pedigree = common::mock_founder_pedigree("parent");
        let ind = pedigree.get_ind_from_label_mut("parent").unwrap();
        ind.set_alleles([0,1]);
        assert!(ind.alleles.is_some());
        ind.clear_alleles();
        assert!(ind.alleles.is_none());


    }

    #[test]

    fn alleles_assignment_founder() {
        let mut pedigree = common::mock_founder_pedigree("parent");
        let iid = pedigree.get_ind_from_label_mut("parent").unwrap().id;
        let result = pedigree.assign_alleles(iid, 0.0, 0, false, &mut Rng::new());
        assert!(result.is_err());
    }

    #[test]
    fn alleles_assignments_unnassigned_parent_alleles(){
        let mut pedigree = common::mock_offspring_pedigree("offspring", None);
        let iid = pedigree.get_ind_from_label_mut("offspring").unwrap().id;
        let result = pedigree.assign_alleles(iid, 0.0, 0, false, &mut Rng::new());
        assert!(result.is_err());
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
        let mut pedigree = common::mock_offspring_pedigree("offspring", None);
        let offspring = pedigree.get_ind_from_label_mut("offspring").unwrap();

        offspring.strands = Some([0,0]);
        let parents_alleles = [[0,1], [0,1]];
        let recombination_prob = 1.0;
        assert_eq!(offspring.currently_recombining, [false, false]);
        let offspring_id = offspring.id;
        perform_allele_assignment(&mut pedigree, offspring_id, parents_alleles, recombination_prob)?;
        let offspring_recombining = pedigree.get_ind(offspring_id).unwrap().currently_recombining;
        assert_eq!(offspring_recombining, [true, true]);
        Ok(())
    }

    //#[test]
    //fn ind_equality() {
    //    let  ind1 = common::mock_founder_pedigree("parent").get_ind_from_label("parent").unwrap().clone();
    //    let  ind2 = common::mock_founder_pedigree("parent").get_ind_from_label("parent").unwrap().clone();
    //    assert_eq!(ind1, ind2);
    //}

    //#[test]
    //fn ind_inequality() {
    //    let  ind1 = common::mock_founder("ind1");
    //    let  ind2 = common::mock_founder("ind2");
    //    assert_ne!(ind1, ind2);
    //}

    //#[test]
    //fn hashable() {
    //    // We're ok here, given that The Hash implementation of Individual only uses the `label` field
    //    // which is not mutable.
    //    #[allow(clippy::mutable_key_type)]
    //    let mut ind_set = HashSet::new();
    //    let n_iters: u32 = 10_000;
    //    for i in 0..n_iters {
    //        let  new_ind = common::mock_founder(&i.to_string());
    //        assert!(ind_set.insert(new_ind.clone()));
    //        assert!(ind_set.contains(&new_ind));
    //    }
    //}

    //#[test]
    //fn ordering() {
    //    let  ind_a = common::mock_founder("A");
    //    let  ind_b = common::mock_founder("B");
    //    assert!(ind_a <  ind_b);
    //    assert!(ind_b <= ind_b);
//
    //    assert!(ind_b >= ind_a);
    //    assert!(ind_b >  ind_a);
//
    //    // Different Index should not impact ordering. What matters is the ID.
    //    let mut ind_a_prime = ind_a.clone();
    //    ind_a_prime.set_tag(SampleTag::new("A", Some(996), None));
    //    assert!(ind_a <= ind_a_prime);
//
//
    //}
//
    //#[test]
    //fn display() {
    //    let (offspring_label, father_label, mother_label) = ("ind1", "ind2", "ind3");
    //    let mut pedigree = common::mock_offspring_pedigree("offspring", Some([father_label, mother_label]));
    //    let offspring = pedigree.get_ind_from_label_mut("offspring").unwrap();
//
    //    let display = format!("{offspring}");
    //    assert!(display.contains(offspring_label));
    //    assert!(display.contains(father_label));
    //    assert!(display.contains(mother_label));
//
    //    offspring.set_tag(SampleTag::new("HG00096", None, None));
    //    let display = format!("{offspring}");
    //    assert!(display.contains("HG00096"));
    //}
}