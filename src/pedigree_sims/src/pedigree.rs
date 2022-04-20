use std::{collections::BTreeMap, path::PathBuf, error::Error};
use crate::io::{Sample, VCFReader};

use std::{
    hash::{Hash, Hasher},
    borrow::Borrow,
    cell::{RefCell, Ref, RefMut},
    rc::Rc,
    cmp::{Ord, Ordering, PartialOrd},
};

#[derive(Debug, Clone)]
pub struct Genome {
    _haplotype1: Vec<u8>,
    _haplotype2: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct Individual {
    pub id   : Option<String>,   //NA005842
    pub label    : String,   //son, mother, stepmom, etc...
    parents  : Option<Parents>,
    _genome   : Genome,
}

type ParentsRef<'a> = (&'a Rc<RefCell<Individual>>, &'a Rc<RefCell<Individual>>);

impl Individual {
    pub fn new(label: String, parents: Option<ParentsRef>) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {id: None, label, parents, _genome: Genome{_haplotype1:vec![], _haplotype2:vec![]}}
    }

    pub fn assign_genome(&mut self, sample: &Sample, vcfs: &Vec<PathBuf> )-> Result<(), Box<dyn Error>> {
        println!("Sample: {}, Id: {}", sample.id(), sample.idx());
        for vcf in vcfs {
            let mut vcf_reader = VCFReader::new(vcf.as_path())?;
            let genome = vcf_reader.parse_sample(sample.idx())?;
            println!("{:?}", genome);
        }
        Ok(())
    }

    pub fn get_parents(&self) -> Option<(Ref<Individual>, Ref<Individual>)> {
        self.parents.as_ref().map(|parents| (RefCell::borrow(&parents.0.0), RefCell::borrow(&parents.0.0)))
    }
    pub fn set_parents(&mut self, parents: ParentsRef) {
        self.parents = Some(Self::format_parents(parents));
    }

    pub fn is_founder(&self) -> bool {
        self.parents == None
    }

    pub fn set_id(&mut self, id: &str){
        self.id = Some(String::from(id));
    }
    fn format_parents(parents:  (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> Parents {
        Parents::new((Rc::clone(parents.0), Rc::clone(parents.1)))
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Parents((Rc<RefCell< Individual>>, Rc<RefCell<Individual>>));

impl Parents{
    pub fn new(parents: (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>)) -> Parents {
        Parents(parents)
    }
}

impl std::fmt::Display for Parents {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} <-> {}", 
            RefCell::borrow(&self.0.0).label,
            RefCell::borrow(&self.0.1).label
        )
    }
}

impl std::fmt::Display for Individual {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let parents = match &self.parents {
            None => "None".to_string(),
            Some(parents) => format!("{}", parents)
        };
        write!(f, "label: {: <10} - parents: {: <25}", self.label, parents)
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

impl Borrow<String> for Individual {
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

#[derive(Debug)]
pub struct Comparison {
    _label           : String,
    pair            : (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>),
    _self_comparison : bool,
}

impl Comparison {
    pub fn new(label: String, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>), self_comparison: bool) -> Comparison {
        let pair = Self::format_pair(pair);
        Comparison{_label: label, pair, _self_comparison: self_comparison }
    }

    pub fn set_pair(&mut self, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)){
        self.pair = Self::format_pair(pair);
    }

    fn format_pair(pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>) {
        (Rc::clone(pair.0), Rc::clone(pair.1))
    }
}

//impl std::fmt::Display for Comparison {
//    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//        unsafe{
//            write!(f, "label: {: <5}\n\t- ind1: ({:?})- ind2: ({:?})", self.label, self.pair.0.borrow(), self.pair.1.borrow())
//        }
//    }
//}

#[derive(Debug)]
pub struct Pedigree {
    pub individuals: BTreeMap<String, Rc<RefCell<Individual>>>,
    comparisons: Vec<Comparison>,
}

impl Pedigree {
    pub fn new() -> Pedigree {
        Pedigree { individuals: BTreeMap::new(), comparisons: Vec::new()}
    }

    pub fn add_individual(&mut self, label: &str, parents: Option<(&String, &String)>) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parents = match parents {
            None => None,
            Some((parent1, parent2)) => {
                let parent1 = self.individuals.get(parent1).ok_or(InvalidInput)?;
                let parent2 = self.individuals.get(parent2).ok_or(InvalidInput)?;
                Some((parent1, parent2))
            },
        };
        let ind = Rc::new(RefCell::new(Individual::new(label.to_owned(), parents)));
        self.individuals.insert(label.to_owned(), ind);
        Ok(())
    }

    pub fn add_comparison(&mut self, label: &str, pair: (&String, &String)) -> std::io::Result<()> {
        use std::io::ErrorKind::InvalidInput;
        let pair0 = self.individuals.get(pair.0).ok_or(InvalidInput)?;
        let pair1 = self.individuals.get(pair.1).ok_or(InvalidInput)?;
        self.comparisons.push(Comparison::new(label.to_owned(), (pair0, pair1), pair0==pair1));
        Ok(())
    }

    pub fn set_relationship(&mut self, ind: &String, parents: (&String, &String)) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parent0 = &self.individuals.get(parents.0).ok_or(InvalidInput)?.clone();
        let parent1 = &self.individuals.get(parents.1).ok_or(InvalidInput)?.clone();

        self.individuals.get_mut(ind)
            .ok_or(InvalidInput)?
            .borrow_mut()
            .set_parents((parent0, parent1));
            Ok(())

    }

    pub fn get_mutind(&mut self, label: &String) -> Result<RefMut<Individual>, std::io::Error>{
        use std::io::ErrorKind::InvalidInput;
        Ok(self.individuals.get_mut(label)
            .ok_or(InvalidInput)?
            .borrow_mut())
    }

    pub fn get_ind(&mut self, label: &String) -> Result<Ref<Individual>, std::io::Error>{
        use std::io::ErrorKind::InvalidInput;
        Ok(RefCell::borrow(self.individuals.get(label).ok_or(InvalidInput)?))
    }


    /// TODO: founders_mut and offsprings_mut have the same code, except the predicate:
    ///  - founters_mut    : ind.is_founder()
    ///  - offsprings_mut : !ind.is_founder()
    pub fn founders_mut(&mut self) -> Vec<RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
            .collect()
    }

    pub fn offsprings_mut(&mut self) -> Vec<RefMut<Individual>> {
        self.individuals
            .values_mut()
            .filter(|ind| !RefCell::borrow(ind).is_founder())
            .map(|x| x.borrow_mut())
            .collect()
    }
    
}

impl Default for Pedigree {
    fn default() -> Self {
        Self::new()
    }
}