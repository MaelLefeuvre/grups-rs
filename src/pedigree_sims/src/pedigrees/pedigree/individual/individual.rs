use std::{
    error::Error,
    hash::{Hash, Hasher},
    cell::RefCell,
    rc::Rc,
    cmp::{Ord, Ordering, PartialOrd},
};

use super::Parents;
use crate::io::vcf::SampleTag;

use genome::Genome;

use rand::{Rng};
use log::{trace};


#[derive(Debug, Clone)]
pub struct Individual {
    pub tag    : Option<SampleTag>,   // NA005842
    pub label  : String,              // son, mother, stepmom, etc...
    parents    : Option<Parents>,
    pub genome : Genome,
    pub strands    : Option<[usize; 2]>,
    pub currently_recombining: [bool; 2],
    pub alleles: Option<[u8; 2]>,
}

const TAG_DISPLAY_LEN: usize = 10;
const LABEL_DISPLAY_LEN: usize = 10;
const PARENTS_DISPLAY_LEN: usize = 25;

impl std::fmt::Display for Individual {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let parents = match &self.parents {
            None => "None".to_string(),
            Some(parents) => format!("{}", parents)
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
    pub fn new(label: &str, parents: Option<ParentsRef>, genome: Genome) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label: label.to_string(), parents, genome, strands: None, currently_recombining: [false, false], alleles: None}
    }

    #[cfg(test)]
    pub fn set_alleles(&mut self, alleles: [u8; 2]) {
        self.alleles = Some(alleles);
    }

    pub fn get_alleles(&self) -> Result<[u8; 2], Box<dyn Error>> {
        self.alleles.ok_or("Attempting to access empty alleles.".into())
    }

    pub fn meiosis(&self, selected_strand: usize, offspring_currently_recombining: bool) -> u8 {
        let selected_strand = match offspring_currently_recombining {
            false => selected_strand,
            true => (selected_strand + 1) % 2
        };
        match self.alleles {
            None          => panic!("Trying to perform meiosis within an empty genome!"),
            Some(alleles) => alleles[selected_strand],
        }
    }

    pub fn assign_strands(&mut self) -> Result<bool, String> {
        if self.parents == None {
            return Err("Assigning strands is meaningless, as this individual has no parents.".to_owned())
        }
        if self.strands != None {
            return Ok(false)
        }
        let mut rng = rand::thread_rng();
        self.strands = Some([rng.gen_range(0, 2), rng.gen_range(0, 2)]);
        Ok(true)
    }

    pub fn get_tag(&self) -> Option<&SampleTag> {
        self.tag.as_ref()
    }

    pub fn set_parents(&mut self, parents: ParentsRef) {
        self.parents = Some(Self::format_parents(parents));
    }

    pub fn is_founder(&self) -> bool {
        self.parents == None
    }

    pub fn set_tag(&mut self, tag: SampleTag){
        self.tag = Some(tag);
    }

    fn format_parents(parents:  [&Rc<RefCell<Individual>>; 2]) -> Parents {
        Parents::new([Rc::clone(parents[0]), Rc::clone(parents[1])])
    }

    pub fn clear_alleles(&mut self){
        self.alleles = None
    }
    
    pub fn assign_alleles(&mut self, recombination_prob: f64, ped_idx: usize) -> Result<bool, Box<dyn Error>> {
        if self.alleles != None {
            return Ok(false)
        }
        match &self.parents {
            None => panic!("Cannot generate genome, as parents are missing."),
            Some(parents) => {
                let mut rng = rand::thread_rng();
                for (i, parent) in parents.iter().enumerate() {

                    //Assign parent genome if not previously generated.
                    if parent.borrow().alleles == None {
                        parent.borrow_mut().assign_alleles(recombination_prob, i).unwrap();
                    }

                    // Check if recombination occured for each parent and update counters if so.
                    if rng.gen::<f64>() < recombination_prob {
                        trace!("Cross-over occured in ped: {:<5} - ind: {}", ped_idx, self.label);
                        self.currently_recombining[i] = ! self.currently_recombining[i];
                    }
                }

                // Assign alleles.
                self.alleles = match self.strands {
                    None => return Err(format!("Cannot assign alleles when self.strands is empty. [{}]", self).into()),
                    Some([s1, s2]) => {
                        let haplo_0 = parents[0].borrow_mut().meiosis(s1, self.currently_recombining[0]);
                        let haplo_1 = parents[1].borrow_mut().meiosis(s2, self.currently_recombining[1]);
                        Some([haplo_0, haplo_1])
                    }
                };
                
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mock_founder(label: &str) -> Individual {
        let genome = Genome::default();
        Individual::new(label, None, genome)
    } 

    fn mock_offspring(label: &str, parents_labels: Option<[&str; 2]>) -> Individual {
        let parents_labels = parents_labels.unwrap_or(["father", "mother"]);

        let genome = Genome::default();
        let father = Rc::new(RefCell::new(mock_founder(parents_labels[0])));
        let mother = Rc::new(RefCell::new(mock_founder(parents_labels[1])));
        let parents = Some([&father, &mother]);
        Individual::new(label, parents, genome)
    }


    fn perform_allele_asignment(offspring: &mut Individual, parents_alleles: [[u8;2];2], recombination_prob: f64) {
        let parents = offspring.parents.as_ref().unwrap();
        parents[0].borrow_mut().alleles = Some(parents_alleles[0]);
        parents[1].borrow_mut().alleles = Some(parents_alleles[1]);        
        offspring.assign_alleles(recombination_prob, 0).unwrap();
    }

    fn run_all_allele_assignment_cases(recombination_prob: f64) {
        let mut offspring = mock_offspring("offspring", None);
        let valid_alleles = vec![[0,0], [0,1], [1,0], [1,1]];
        let mut valid_strands = vec![[0,0], [0,1], [1,0], [1,1]];

        for parent_0_alleles in valid_alleles.iter(){
            for parent_1_alleles in valid_alleles.iter() {
                for strands in valid_strands.iter_mut() {
                    offspring.strands = Some(*strands);
                    perform_allele_asignment(&mut offspring, [*parent_0_alleles, *parent_1_alleles], recombination_prob);

                    // If the individual's parent is 'recombining', we expect strand assignment to be inverted. 0 becomes 1 ; 1 becomes 0
                    for i in [0,1]{
                        if offspring.currently_recombining[i] {
                            strands[i] = (strands[i]+1) % 2;
                        }
                    }

                    let want = [parent_0_alleles[strands[0]], parent_1_alleles[strands[1]]];
                    assert_eq!(offspring.alleles.unwrap(), want);
                    offspring.alleles = None;
                }
            }
        }
    }

     #[test]
    fn alleles_getter_filled(){
        let mut ind = mock_founder("offspring");
        let alt_ref = [0,1];
        for i in alt_ref {
            for j in alt_ref {
               let alleles = [i, j];
               ind.set_alleles(alleles);
               assert_eq!(ind.get_alleles().unwrap(), alleles);
            }
        }
    }

    #[test]
    fn alleles_getter_empty(){
        let ind = mock_founder("offspring");
        let alleles = ind.get_alleles();
        assert!(alleles.is_err());
    }

    #[test]
    fn meiosis_not_recombining() {
        let offspring_currently_recombining = false; 
        let mut ind = mock_founder("offspring");
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
        let mut ind = mock_founder("parent");
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
        let ind = mock_founder("parent");
        ind.meiosis(0, false);
    }


    #[test]
    fn strand_setter_offspring() {
        let valid_strands = vec![[0,0], [0,1], [1,0], [1,1]];
        let mut ind = mock_offspring("offspring", None);

        for _ in 0..1000 {
            ind.assign_strands().unwrap();
            let rng_strands = ind.strands.unwrap();
            assert!(valid_strands.contains(&rng_strands));
        }
    }

    #[test]
    fn strand_setter_founder() {
        let mut ind = mock_founder("parent");
        let result = ind.assign_strands();
        assert!(result.is_err())
    }

    #[test]
    fn get_set_sampletag() {
        let mut ind = mock_founder("parent");
        let tag = SampleTag::new("HG00096", Some(0));
        ind.set_tag(tag.to_owned());
        assert_eq!(ind.get_tag(), Some(&tag));
        
    }

    #[test]
    fn get_empty_tag(){
       let ind = mock_founder("parent");
       let result = ind.get_tag();
       assert!(result.is_none());
    }

    #[test]
    fn founder_is_founder() {
        let ind = mock_founder("parent");
        assert!(ind.is_founder());
    }

    #[test]
    fn offpsring_is_not_founder() {
        let ind = mock_offspring("offspring", None);
        assert!(!ind.is_founder());
    }

    #[test]
    fn clear_alleles() {
        let mut ind = mock_founder("parent");
        ind.set_alleles([0,1]);
        assert!(ind.alleles.is_some());
        ind.clear_alleles();
        assert!(ind.alleles.is_none());


    }

    #[test]
    #[should_panic]
    fn alleles_assignment_founder() {
        let mut ind = mock_founder("parent");
        let _result = ind.assign_alleles(0.0, 0);
        //assert!(result.is_err());
    }

    #[test]
    #[should_panic]
    fn alleles_assignments_unnassigned_parent_alleles(){
        let mut ind = mock_offspring("offspring", None);
        let _result = ind.assign_alleles(0.0, 0);
    }

    #[test]
    fn alleles_assignment_offspring_no_recombination() {
        run_all_allele_assignment_cases(0.0);
    }
    #[test]
    fn alleles_assignment_offspring_full_recombination() {
        run_all_allele_assignment_cases(1.0);
    }

    #[test]
    fn allele_assignment_updates_recombination_status(){
        let mut offspring = mock_offspring("offspring", None);
        offspring.strands = Some([0,0]);
        let parents_alleles = [[0,1], [0,1]];
        let recombination_prob = 1.0;

        assert_eq!(offspring.currently_recombining, [false, false]);
        perform_allele_asignment(&mut offspring, parents_alleles, recombination_prob);
        assert_eq!(offspring.currently_recombining, [true, true]);


    }

    #[test]
    fn ind_equality() {
        let  ind1 = mock_founder("parent");
        let  ind2 = mock_founder("parent");
        assert_eq!(ind1, ind2)
    }

    #[test]
    fn ind_inequality() {
        let  ind1 = mock_founder("ind1");
        let  ind2 = mock_founder("ind2");
        assert_ne!(ind1, ind2)
    }

    #[test]
    fn hashable() {
        let mut ind_set = std::collections::HashSet::new();
        let n_iters: u32 = 10_000;
        for i in 0..n_iters {
            let  new_ind = mock_founder(&i.to_string());
            assert!(ind_set.insert(new_ind.clone()));
            assert!(ind_set.get(&new_ind).is_some());
        }
    }

    #[test]
    fn ordering() {
        let  ind_a = mock_founder("A");
        let  ind_b = mock_founder("B");
        assert!(ind_a <  ind_b);
        assert!(ind_b <= ind_b);

        assert!(ind_b >= ind_a);
        assert!(ind_b >  ind_a);

        // Different Index should not impact ordering. What matters is the ID.
        let mut ind_a_prime = ind_a.clone();
        ind_a_prime.set_tag(SampleTag::new("A", Some(996)));
        assert!(! (ind_a >  ind_a_prime));


    }

    #[test]
    fn display() {
        let (offspring_label, father_label, mother_label) = ("ind1", "ind2", "ind3");
        let mut offspring = mock_offspring(offspring_label, Some([father_label, mother_label]));
        
        let display = format!("{offspring}");
        assert!(display.contains(offspring_label));
        assert!(display.contains(father_label));
        assert!(display.contains(mother_label));

        offspring.set_tag(SampleTag::new("HG00096", None));
        let display = format!("{offspring}");
        assert!(display.contains("HG00096"));
    }
}