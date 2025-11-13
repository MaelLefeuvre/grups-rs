//mod mod_copy;

use std::fmt::{self, Display, Formatter};

use crate::pedigrees::pedigree::individual::IndividualError;
use crate::pedigrees::constants::{AVG_PWD_FORMAT_LEN, COMPARISON_LABEL_FORMAT_LEN, FLOAT_FORMAT_PRECISION, IND_LABEL_FORMAT_LEN, IND_TAG_FORMAT_LEN, OVERLAP_FORMAT_LEN, PWD_FORMAT_LEN, SEX_FORMAT_LEN};

mod relationship;
use fastrand::Rng;
use itertools::Itertools;
use log::trace;
use relationship::{Relationship, RelationshipId};

mod individual;
use individual::{Individual, IndividualId};

mod comparisons;
use comparisons::{PedComparisons, PedComparison};


use slotmap::SlotMap;
use std::collections::BTreeMap;

use located_error::prelude::*;
use grups_io::read::{
    genotype_reader::GenotypeReader, PanelReader, SampleTag
};

use genome::Sex;

use fastrand;


#[cfg(test)] mod tests ; 

pub mod parser;

mod contaminant;
pub use contaminant::Contaminant;



pub mod pedparam;
use pedparam::PedigreeParams;

mod error;
pub use error::PedigreeError;


// --------------------------------------------------------------------------------------------- //
// Individuals

#[derive(Debug, Clone)]
pub struct PedIndividuals {
    inner: SlotMap<IndividualId, Individual>,
    labels: BTreeMap<String, IndividualId>,
}

impl Default for PedIndividuals {
    fn default() -> Self { Self::new() }
}

impl PedIndividuals {
    pub fn new() -> Self {
        Self{inner: SlotMap::with_key(), labels: BTreeMap::new()}
    }

    pub fn get_ind_id(&self, label: &str) -> Option<IndividualId> {
        self.labels.get(label).copied()
    }

    pub fn get_ind(&self, id: IndividualId) -> Option<&Individual> {
        self.inner.get(id)
    }

    pub fn get_ind_from_label_mut(&mut self, label: &str) -> Option<&mut Individual> {
        let id = self.labels.get(label)?;
        self.get_ind_mut(*id)
    }

    pub fn get_ind_from_label(&self, label: &str) -> Option<&Individual> {
        let id = self.labels.get(label)?;
        self.get_ind(*id)
    }

    pub fn get_ind_mut(&mut self, id: IndividualId) -> Option<&mut Individual> {
        self.inner.get_mut(id)
    }


    // ---- NB: Sorting is only required to ensure backwards compatibility with the current test-suite. 
    // TODO: remove this alltogether.
    #[inline] pub fn _sorted_iter<P: FnMut(&&Individual) -> bool>(&self, predicate: P) -> impl Iterator<Item = &Individual> {
        Itertools::sorted_by(self.inner.values().filter(predicate), |a,b| a.label().partial_cmp(b.label()).unwrap())
    }
    
    // ---- NB: Sorting is only required to ensure backwards compatibility with the current test-suite. 
    // TODO: remove this alltogether.
    #[inline] pub fn _sorted_iter_mut<T: FnMut(&&mut Individual) -> bool>(&mut self, predicate: T) -> impl Iterator<Item = &mut Individual> {
        Itertools::sorted_by(self.inner.values_mut().filter(predicate), |a,b| a.label().partial_cmp(b.label()).unwrap())
    }

    #[inline]
    pub fn founders_mut(&mut self) -> impl Iterator<Item = &mut Individual> {
        self.inner.values_mut().filter(|ind| ind.is_founder())      // Unsorted
    }
    
    #[inline]
    pub fn offsprings_ids(&self) -> Vec<IndividualId> {
        self.inner.values().filter(|ind| ind.is_offspring()).map(|ind| ind.id).collect::<Vec<IndividualId>>()  
    }

    // NB: This is only used to ensure backwards compatibility with the current test-suite, and should be replaced ASAP.
    #[inline]
    pub fn offsprings_ids_sorted(&self) -> Vec<IndividualId> { 
        self._sorted_iter(|ind| ind.is_offspring()).map(|ind| ind.id).collect::<Vec<IndividualId>>()    // Sorted (inefficient)
    }

    #[inline]
    pub fn offsprings(&self) -> impl Iterator<Item = &Individual> {
        self.inner.values().filter(|ind| ind.is_offspring())           // Unsorted
    }

    #[inline]
    pub fn offsprings_mut(&mut self) -> impl Iterator<Item = &mut Individual> {
        self.inner.values_mut().filter(|ind| ind.is_offspring())      // Unsorted
    }

    pub fn add_individual(&mut self, individual: &str, parents: Option<[&str; 2]>, sex:Option<Sex>) -> IndividualId {
        // ---- Add parents as individuals if there were not previously found
        let _parent_ids = parents.map(|parents| {
            parents.map(|parent| {
                self.get_ind_id(parent)
                    .or_else(|| Some(self.add_individual(parent, None, None))) 
                    .expect("Individual should be includable")
            })
        });
        // ---- Insert individuals, as well as its parents
        // ---- Early return if ind already exists.
        let ind_id = if let Some(id) = self.labels.get(individual) { *id } else {
            let id = self
                .inner
                .insert_with_key(|id| Individual::new(id, individual, sex));
            self.labels.insert(individual.to_string(), id);
            id
        };
        ind_id
    }

    /// Set all individual allele's to `None` within this pedigree.
    pub fn clear_alleles(&mut self){
        self.inner.values_mut().for_each(|ind| ind.clear_alleles())
    }
}



#[derive(Debug, Clone)]
pub struct Pedigree {
    pub individuals: PedIndividuals,
    pub edges: SlotMap<RelationshipId, Relationship>,
    pub comparisons: PedComparisons,
    params: Option<PedigreeParams>,
    pop: Option<String>,
}

impl Default for Pedigree {
    fn default() -> Self { Self::new() }
}

impl Pedigree {
    pub fn new() -> Self {
        Self {
            individuals: PedIndividuals::new(),
            edges: SlotMap::with_key(),
            comparisons: PedComparisons::new(),
            params: None,
            pop: None,
        }
    }

    #[inline]
    pub fn clear_alleles(&mut self) {
        self.individuals.clear_alleles();
    }

    #[inline]
    pub fn compare_alleles(&mut self, contam_pop_af: [f64; 2], pileup_error_probs: &[f64; 2], rng: &mut Rng) -> Result<()> {
        use PedigreeError::FailedAlleleComparison;
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
            let alleles = comparison.pair.map(|ind_id| {
                self.individuals.get_ind(ind_id)
                    .expect("Individual should exist at this point")
                    .alleles
                    .expect("Alleles should be set at this point")
            });
            comparison.compare_alleles(alleles, contam_rate, contam_pop_af, seq_error_rate, rng).with_loc(|| FailedAlleleComparison)?;
        }
        Ok(())
    }

    #[inline]
    pub fn update_founder_alleles(&mut self, reader: &dyn GenotypeReader, rng: &mut Rng) -> Result<()> {
        let loc_msg = "While updating founder individiuals' alleles";
        // ---- Extract this pedigree allele frequency downsampling rate.
        let af_downsampling_rate = self.get_params().loc(loc_msg)?.af_downsampling_rate;

        // ---- Perform allele fixation at random, according to this pedigrees af_downsampling_rate
        if rng.f64() < af_downsampling_rate {
            self.individuals.founders_mut().for_each(|founder| founder.alleles = Some([0, 0]));
        } else {
            // ---- Fetch and assign the 'true' alleles for all founder individuals. 
            for founder in self.individuals.founders_mut() {
                // ---- Extract founder tag ; raise an error if None is returned.
                let founder_tag = founder.get_tag().loc(loc_msg)?;

                // ---- Fetch and assign the individual's allele from our reader.
                founder.alleles = Some(reader.get_alleles(founder_tag)?);
            }
        }
        Ok(())
    }

    #[inline]
    pub fn compute_offspring_alleles(&mut self, interval_prob_recomb: f64, pedigree_index: usize, xchr_mode: bool, rng: &mut fastrand::Rng) -> Result<()> {
        for offspring_id in self.individuals.offsprings_ids_sorted() {
            self.assign_alleles(offspring_id, interval_prob_recomb, pedigree_index, xchr_mode, rng)
                .with_loc(|| format!("While attempting to assign the alleles of {}", 
                    self.individuals.get_ind(offspring_id).expect("Individual should be retrievable").label())
                )?;
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
    
    pub fn assign_offspring_strands(&mut self) -> Result<()> {
        for offspring_id in self.individuals.offsprings_ids_sorted() {
            let ind = self.individuals.get_ind_mut(offspring_id).expect("Individual should not be retrievable.");
            ind.assign_strands().with_loc(||PedigreeError::FailedAlleleAssignment(ind.label().to_string()))?;
        }
        Ok(())
    }

    /// Randomly assign the sex of each individual.
    pub fn assign_random_sexes(&mut self) -> Result<()> {
        for offspring_id in self.individuals.offsprings_ids_sorted() {
            self.assign_random_sex(offspring_id).with_loc(||PedigreeError::FailedSexAssignment(
                self.individuals.get_ind(offspring_id).expect("Individual should be retrievable").label().to_string()
            ))?;
        }
        Ok(())
    }

    pub fn add_individual(&mut self, individual: &str, parents: Option<[&str; 2]>, sex:Option<Sex>) -> IndividualId {
        // ---- Add parents as individuals if there were not previously found
        let parent_ids = parents.map(|parents| {
            parents.map(|parent| {
                self.individuals.get_ind_id(parent)
                    .or_else(|| Some(self.add_individual(parent, None, None))) 
                    .expect("Individual should be includable")
            })
        });

        let ind_id = self.individuals.add_individual(individual, parents, sex);

        if let Some(parents) = parent_ids {
            if !self.individuals.get_ind(ind_id).expect("Individual should be retrievable").has_parents() {
                // ---- Create relationship from ind to parent 1 and parent 2
                println!("----- While adding individual ! ({individual})");
                self.set_relationship(ind_id, parents);
            }
        }
        ind_id
    }

    pub fn add_comparison(
        &mut self,
        label: &str,
        pair: [&str; 2],
    ) -> std::io::Result<()> {
        use std::io::ErrorKind::InvalidInput;
        // @TODO: converting to a Vec<T> and back to a [T; 2] is not that elegant, but
        // [].try_map() is still unstable...
        let pair_ids: [IndividualId; 2] = pair
            .iter()
            .flat_map(|ind_label| self.individuals.get_ind_id(ind_label).ok_or(InvalidInput))
            .collect::<Vec<_>>()
            .try_into()
            .expect("Vec should be of length 2");
        self.comparisons.push(PedComparison::new(label, pair_ids));
        Ok(())
    }

    pub fn set_relationship(
        &mut self,
        from: IndividualId,
        to: [IndividualId; 2],
    ) -> [RelationshipId; 2] {
        //println!("  - Setting relationship: {from} -> {to:?}");
        let parental_relationships: [RelationshipId; 2] = to.map(|parent_id| {
            self.edges
                .insert_with_key(|id| Relationship::new(id, from, parent_id))
        });

        if let Some(ind) = self.individuals.get_ind_mut(from) {
            ind.add_parents(parental_relationships);
        }

        parental_relationships
    }

    pub fn set_founder_tags(&mut self, panel: &PanelReader, pop: &String, contaminants: Option<&Contaminant>) -> Result<()>{
        use PedigreeError::{MissingContaminant, MissingSampleTag};
        self.pop = Some(pop.to_owned());

        // ---- Contaminating individual are excluded from pedigree individuals.
        let mut exclude_tags = contaminants.map(|cont| cont.as_flat_list()).with_loc(||MissingContaminant)?;

        // ---- For each founder, pick and assign a random SampleTag using our panel (without replacement)
        for founder in self.individuals.founders_mut() {
            // ---- Pick a random SampleTag from our panel.
            let random_tag = panel.random_sample(pop, Some(&exclude_tags), founder.sex)?.with_loc(||MissingSampleTag)?; 
            founder.set_tag(random_tag.clone());
            exclude_tags.push(random_tag); // Exclude this pedigree individual for other iterations.
        };
        Ok(())
    }

    pub fn get_parents_ids(&self, id: IndividualId) -> Option<[IndividualId; 2]> {
        if let Some(ind) = self.individuals.get_ind(id) {
            if let Some(parent_ids) = ind.get_parents() {
                let parents_ids = parent_ids.map(|rel_id| {
                    self.edges.get(rel_id).expect("Relationship should be retrievable").to
                });
                return Some(parents_ids)
            }
        }
        None
    }

    pub fn get_parents(&self, id: IndividualId) -> Option<[&Individual; 2]> {
        self.get_parents_ids(id).map(|ids| {
            ids.map(|parent_id| self.individuals.get_ind(parent_id).expect("Individual should be retrievable"))
        })
    }


    #[inline]
    pub fn assign_alleles (&mut self, iid: IndividualId, recombination_prob: f64, ped_idx: usize, xchr_mode: bool, rng: &mut Rng) -> Result<bool> {
        use IndividualError::{InvalidAlleleAssignment, MissingParents, MissingStrands};
        // ---- Ensure this method call is non-redundant.
        if self.individuals.get_ind(iid).expect("Individual should be retrievable").alleles.is_some() {
            return Ok(false)
        }

        // ---- Ensure the individual is 'equipped' with parents and strands before attempting allele assignment.
        let Some(parents) = &self.get_parents_ids(iid) else {
            return Err(anyhow!(MissingParents)).with_loc(||InvalidAlleleAssignment)
        };

        // ---- perform allele assignment for each parent
        for (i, parent_id) in parents.iter().enumerate() {
            // ---- Assign parent genome if not previously generated
            let parent = self.individuals.get_ind(*parent_id).expect("Individual should be retrievable");
            if parent.alleles.is_none() {
                self.assign_alleles(*parent_id, recombination_prob, ped_idx, xchr_mode, rng)?;
            }

            // ---- Check if recombination occured for each parent and update recombination tracker if so
            let parent = self.individuals.get_ind(*parent_id).expect("Individual should be retrievable");
            if (!xchr_mode || parent.sex == Some(Sex::Female)) && rng.f64() < recombination_prob {
                let ind = self.individuals.get_ind(iid).expect("Individual should be retrievable");
                trace!("- Cross-over occured in ped: {:<5} - ind: {} ({} {:?})", ped_idx, ind.label(), parent.label(), parent.sex);
                let ind = self.individuals.get_ind_mut(iid).expect("Individual should be retrievable"); // TODO: find a way to appease the borrow checker and remove this redundant borrow
                ind.currently_recombining[i] = ! ind.currently_recombining[i];
            }
        }

        // ---- Perform allele assignment for `ind`, by simulating meiosis for each parent.
        let ind = self.individuals.get_ind(iid).expect("Individual should be retrievable");
        let Some(strands) = ind.strands else {
            return Err(anyhow!(MissingStrands)).with_loc(||InvalidAlleleAssignment)
        };

        let ind_alleles = if xchr_mode {
            let mut alleles = [0u8; 2];
            // ---- Find the index of both parents
            let father_idx = parents.iter().position(|p| self.individuals.get_ind(*p).expect("Parent should be retrievable").sex == Some(Sex::Male)).expect("No parent found..");
            let mother_idx = (father_idx + 1 ) % 2;

            // -- assign maternal strand
            let mat_strand_currently_recombining = ind.currently_recombining[mother_idx];
            alleles[mother_idx] = self.individuals.get_ind(parents[mother_idx]).expect("Parent should be retrievable").meiosis(strands[mother_idx], mat_strand_currently_recombining);

            // -- Assign paternal strand
            let pat_strand_currently_recombining = ind.currently_recombining[father_idx];
            alleles[father_idx] = match ind.sex {
                Some(Sex::Male)           => Ok(alleles[mother_idx]), // If the descendant is a male, alleles are exclusively from the mother
                Some(Sex::Female)         => Ok(self.individuals.get_ind(parents[father_idx]).expect("Parent should be retrievable").meiosis(strands[father_idx], pat_strand_currently_recombining)),
                Some(Sex::Unknown) | None => Err(IndividualError::UnknownOrMissingSex).loc("While attempting to assign alleles during X-chromosome-mode"),
            }?;

            // ---- Sanity checks
            if ind.sex == Some(Sex::Male) && alleles[0] != alleles[1] {
                // Male individuals are not expected to be heterozygous during X-chromosome simulations.
                return Err(IndividualError::SpuriousAlleleAssignment{alleles})
                    .loc("Male Individual is heterozygous while in X-chromosome mode")
            }

            if ind.currently_recombining[father_idx] {
                // Fathers are not expected to recombine during X-chromosome simulations.
                return Err(IndividualError::InvalidOrSpuriousRecombinationEvent)
                    .loc("Father is recombining while in X-chromosome-mode")
            }

            Some(alleles)
        } else {
            let haplo_0 = self.individuals.get_ind(parents[0]).expect("Parent should be retrievable").meiosis(strands[0], ind.currently_recombining[0]);
            let haplo_1 = self.individuals.get_ind(parents[1]).expect("Parent should be retrievable").meiosis(strands[1], ind.currently_recombining[1]);
            Some([haplo_0, haplo_1])
        };

        self.individuals.get_ind_mut(iid).expect("Individual should be retrievable").alleles = ind_alleles;
        Ok(true)
    }

    #[inline]
    pub fn assign_random_sex(&mut self, id: IndividualId) -> Result<bool> {
        use IndividualError::InvalidSexAssignment;
        // ---- Ensure this method call is non-redundant
        if self.individuals.get_ind(id).expect("Individual should be retrievable").sex.is_some() {
            return Ok(false)
        }

        // ---- Perform sex-asignment for each parent (if the individual has known parents.)
        if let Some(parents_ids) = self.get_parents_ids(id) {
            for (i, parent_id) in parents_ids.iter().enumerate() {
                // ---- Assign parent sex if not previously decided
                let spouse_id  = parents_ids[(i+1) % 2];
                let parent_sex =  self.individuals.get_ind(*parent_id).expect("Individual should be retrievable").sex;
                let spouse_sex =  self.individuals.get_ind(spouse_id).expect("Individual should be retrievable").sex;
                if parent_sex.is_none() { // If the parent's sex is still unknown
                    let parent_sex = match spouse_sex { 
                        Some(sex) => match sex {  // If the spouse sex is already known, assign the opposite sex.
                            Sex::Female => Some(Sex::Male),
                            Sex::Male   => Some(Sex::Female),
                            Sex::Unknown => None
                        }
                        None => Some(Sex::random()), // If not, assign a random sex to the parent.
                    };
                    let parent = self.individuals.get_ind_mut(*parent_id).expect("Individual should be retrievable");
                    parent.sex = parent_sex;
                }
                // ---- Apply the same process for the parent
                self.assign_random_sex(*parent_id).with_loc(|| InvalidSexAssignment)?;
            }
        }

        // ---- Randomly assign sex of the considered individual
        let ind = self.individuals.get_ind_mut(id).expect("Individual should be retrievable");
        ind.sex = Some(Sex::random());

        Ok(true)
    }

    pub fn all_sex_assigned(&self) -> bool {
        self.individuals.inner.values().all(|ind| ind.is_sex_assigned())
    }
    
    pub fn _display_comparison(&self, f: &mut Formatter<'_>, comp: &PedComparison) -> fmt::Result {
            let default_tag = SampleTag::new("None", None, None);
            let ind1 = self.individuals.get_ind(comp.pair[0]).expect("Individual should be retrievable");
            let ind2 = self.individuals.get_ind(comp.pair[1]).expect("Individual should be retrievable");

            // <Comparison-Label> <Ind1-label> <Ind2-label> <Ind1-reference> <Ind2-reference> <Sum.PWD> <Overlap> <Avg.PWD> <Ind1-sex> <Ind2-sex>
            writeln!(f,
                "{: <COMPARISON_LABEL_FORMAT_LEN$} - \
                {: <IND_LABEL_FORMAT_LEN$} - \
                {: <IND_LABEL_FORMAT_LEN$} - \
                {: <IND_TAG_FORMAT_LEN$} - \
                {: <IND_TAG_FORMAT_LEN$} - \
                {: >PWD_FORMAT_LEN$} - \
                {: >OVERLAP_FORMAT_LEN$} - \
                {: <AVG_PWD_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$} - \
                {: <SEX_FORMAT_LEN$} - \
                {: <SEX_FORMAT_LEN$}",
                comp.label,
                ind1.label(),
                ind2.label(),
                ind1.get_tag().as_ref().unwrap_or(&&default_tag).id(),
                ind2.get_tag().as_ref().unwrap_or(&&default_tag).id(),
                comp.get_sum_pwd(),
                comp.get_overlap(),
                comp.get_avg_pwd(),
                ind1.sex.map_or(String::from("None"), |s| s.to_string()),
                ind2.sex.map_or(String::from("None"), |s| s.to_string())
            )
    }
}

impl Display for Pedigree {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        self.comparisons.iter().try_fold((), |_, comp| {
            self._display_comparison(f, comp)
        })
    }
}


#[cfg(test)]
mod test {
    use super::*;

    fn test_pedigree_set() -> Pedigree {
        let mut pedigree = Pedigree::new();
        pedigree.add_individual("father", None, None);
        pedigree.add_individual("mother", None, None);
        pedigree.add_individual("offspr", Some(["father", "mother"]), None);
    
        let father = pedigree.individuals.get_ind_from_label_mut("father").expect("Cannot extract father");
        father.set_alleles([0, 1]);
    
        let mother = pedigree.individuals.get_ind_from_label_mut("mother").expect("Cannot extract mother");
        mother.set_alleles([1, 0]);
    
    
        pedigree
    }

    type TestPedDef = Vec<(&'static str, Option<[&'static str; 2]>)>;
    
    fn test_pedigree_random(def: Option<TestPedDef>) -> std::io::Result<Pedigree> {
        let mut pedigree = Pedigree::new();
        if let Some(map) = def {
            for (label, parents) in map {
                println!("{label}, {parents:?}");
                pedigree.add_individual(label, parents, None);
            }
            for ind in pedigree.individuals.inner.values_mut() {
                if ind.is_founder() {
                    ind.set_alleles([u8::from(fastrand::bool()), u8::from(fastrand::bool())]);
                }
            }
        } else {
            let def = Vec::from([
                ("father", None),
                ("mother", None),
                ("offspr", Some(["father", "mother"]))
            ]);
            return test_pedigree_random(Some(def))
        }

        Ok(pedigree)
    }


    #[test]
    #[should_panic = "Failed to assign alleles"]
    fn meiosis_assign_alleles_empty_strands(){
        let mut pedigree = test_pedigree_set();
        let offspr = pedigree.individuals.get_ind_id("offspr").expect("Cannot extract offspr");
        pedigree.assign_alleles(offspr, 0.0, 0, false, &mut fastrand::Rng::new()).expect("Failed to assign alleles");
    }

    #[test]
    fn meiosis_assign_alleles_filled_strands(){
        let mut rng      = fastrand::Rng::new();
        let mut pedigree = test_pedigree_set();
        let offspr       = pedigree.individuals.get_ind_from_label_mut("offspr").expect("Cannot extract offspr");
        offspr.strands   = Some([0, 0]);
        
        let offspr_id   = pedigree.individuals.get_ind_id("offspr").expect("Cannot extract offspr");
        let output = pedigree.assign_alleles(offspr_id, 0.0, 0, false, &mut rng).expect("Failed to assign alleles");
        assert!(output);
        let output = pedigree.assign_alleles(offspr_id, 0.0, 0, false, &mut rng).expect("Failed to assign alleles");
        assert!(!output);
    }

    #[test]
    fn meiosis_check_strands_00() {
        let mut pedigree = test_pedigree_set();
        let offspr_id    = pedigree.individuals.get_ind_id("offspr").expect("Cannot extract offspr");
        let offspr   = pedigree.individuals.get_ind_mut(offspr_id).expect("Cannot extract offspr");
        offspr.strands   = Some([0, 0]);
        pedigree.assign_alleles(offspr_id, 0.0, 0, false, &mut fastrand::Rng::new()).expect("Failed to assign alleles");
        assert_eq!(pedigree.individuals.get_ind(offspr_id).expect("Individual should be retrievable").alleles, Some([0, 1]));
    }

    #[test]
    fn meiosis_check_strands_01() {
        let mut pedigree = test_pedigree_set();
        let offspr_id    = pedigree.individuals.get_ind_id("offspr").expect("Cannot extract offspr");
        let offspr       = pedigree.individuals.get_ind_mut(offspr_id).expect("Cannot extract offspr");
        offspr.strands   = Some([1, 1]);

        pedigree.assign_alleles(offspr_id, 0.0, 0, false, &mut fastrand::Rng::new()).expect("Failed to assign alleles");
        assert_eq!(pedigree.individuals.get_ind(offspr_id).expect("Individual should be retrievable").alleles, Some([1, 0]));
        assert_eq!(pedigree.individuals.get_ind(offspr_id).expect("Individual should be retrievable").currently_recombining, [false, false]);

    }

    #[test]
    fn meiosis_check_recombination() {
        let mut pedigree = test_pedigree_set();
        let offspr_id    = pedigree.individuals.get_ind_id("offspr").expect("Cannot extract offspr");
        let offspr       = pedigree.individuals.get_ind_mut(offspr_id).expect("Cannot extract offspr");
        offspr.strands   = Some([0, 1]);
        pedigree.assign_alleles(offspr_id, 1.0, 0, false, &mut fastrand::Rng::new()).expect("Failed to assign alleles");

        assert_eq!(pedigree.individuals.get_ind(offspr_id).expect("Individual should be retrievable").alleles, Some([1, 1]));
        assert_eq!(pedigree.individuals.get_ind(offspr_id).expect("Individual should be retrievable").currently_recombining, [true, true]);
    }

    #[test]
    fn assign_offspring_strands() -> Result<()> {
        let mut pedigree = test_pedigree_set();
        pedigree.assign_offspring_strands()?;
        Ok(())
    }

    #[test]
    fn assign_sexes() -> Result<()> {
        let def = Vec::from([
            ("F1.1", None),
            ("F1.2", None),
            ("F1.3", None),
            ("F2.1", None),
            ("O2.2", Some(["F1.1", "F1.2"])),
            ("O2.3", Some(["F1.1", "F1.2"])),
            ("F2.4", None),
            ("02.5", Some(["F1.2", "F1.3"])),
            ("O3.1", Some(["F2.1", "O2.2"])),
            ("O3.2", Some(["O2.3", "F2.4"])),        
        ]);
        for _ in 0..1000 {
            let mut pedigree = test_pedigree_random(Some(def.clone())).expect("Pedigree should be constructible");
            pedigree.assign_random_sexes()?;
            for ind in pedigree.individuals.inner.values() {
                assert!(ind.sex.is_some());
            }
        }
        Ok(())
    }

    #[test]
    fn clear_alleles() {
        let mut pedigree = test_pedigree_set();
        let mut rng = Rng::new();
        for ind in pedigree.individuals.inner.values_mut() {
            ind.alleles = Some([u8::from(rng.bool()), u8::from(rng.bool())])
        }
        assert!(pedigree.individuals.inner.values().all(|ind| ind.alleles.is_some()));
        pedigree.clear_alleles();
        assert!(pedigree.individuals.inner.values().all(|ind| ind.alleles.is_none()));

    }
}
