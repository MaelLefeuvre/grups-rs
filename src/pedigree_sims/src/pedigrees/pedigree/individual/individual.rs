use std::{ error::Error, };


use crate::{
    io::{
        vcf::{
            SampleTag,
        }, 
    }
};

use genome::{
    Genome,
};


use rand::{Rng};
use log::{trace};



use std::{
    hash::{Hash, Hasher},
    cell::{RefCell},
    rc::Rc,
    cmp::{Ord, Ordering, PartialOrd},
};

use super::Parents;

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
        write!(f, "tag: {: <10} label: {: <10} - parents: {: <25}", tag, self.label, parents)
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
    pub fn new(label: String, parents: Option<ParentsRef>, genome: Genome) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label, parents, genome, strands: None, currently_recombining: [false, false], alleles: None}
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
