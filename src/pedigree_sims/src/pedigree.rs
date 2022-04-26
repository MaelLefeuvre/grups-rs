use std::{collections::{BTreeMap, HashSet}, path::PathBuf, error::Error};
use crate::io::{SampleTag, VCFReader, VCFPanelReader};
use log::{info};
use pwd_from_stdin::genome::{SNPCoord, Genome, GeneticMap, Chromosome};
use rand::Rng;
use rand::seq::SliceRandom;



use std::{
    hash::{Hash, Hasher},
    cell::{RefCell, Ref, RefMut},
    rc::Rc,
    cmp::{Ord, Ordering, PartialOrd},
};


#[derive(Debug, Clone)]
pub struct Individual {
    pub tag    : Option<SampleTag>,   // NA005842
    pub label  : String,              // son, mother, stepmom, etc...
    parents    : Option<Parents>,
    pub genome : Genome,
}

type ParentsRef<'a> = (&'a Rc<RefCell<Individual>>, &'a Rc<RefCell<Individual>>);

impl Individual {
    pub fn new(label: String, parents: Option<ParentsRef>, genome: Genome) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label, parents, genome}
    }

    pub fn get_tag(&self) -> Option<&SampleTag> {
        self.tag.as_ref()
    }

    pub fn add_locus(&mut self, chromosome: &u8, pos: u32, alleles: (u8, u8), af: f64) -> Result<(), &str> {
        self.genome.get_chr_mut(chromosome).ok_or("Missing chromosome.")?.add_locus(pos, alleles, af);
        Ok(())
    }

    pub fn get_parents(&self) -> Option<(Ref<Individual>, Ref<Individual>)> {
        self.parents.as_ref().map(|parents| (RefCell::borrow(&parents.0), RefCell::borrow(&parents.0)))
    }

    pub fn get_parents_mut(&self) -> Option<(RefMut<Individual>, RefMut<Individual>)> {
        self.parents.as_ref().map(|parents| (RefCell::borrow_mut(&parents.0), RefCell::borrow_mut(&parents.0)))
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

    fn format_parents(parents:  (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> Parents {
        Parents::new((Rc::clone(parents.0), Rc::clone(parents.1)))
    }

    pub fn has_empty_genome(&self) -> bool {
        self.genome.is_empty()
    }

    pub fn generate_genome(&mut self, genetic_map: &GeneticMap) -> Result<bool, Box<dyn Error>>{
        if ! self.has_empty_genome() {
            println!("{} Genome already generated.", self.label);
            return Ok(false)
        }

        println!("Generating {} genome.", self.label);
        match &self.parents {
            None => panic!("Cannot generate genome, as parents are missing."),
            Some(parents) => {
                for parent in [&parents.0, &parents.1] {
                    if parent.borrow().has_empty_genome() {
                        println!(" -> Parent {} has empty genome. Generating...", parent.borrow().label);
                        parent.borrow_mut().generate_genome(genetic_map).unwrap();
                    }
                    else {
                        println!("-> Parent {} has genome.", parent.borrow().label);
                    }
                }
                let gamete_1 = parents.0.borrow_mut().genome.meiosis(genetic_map);
                let gamete_2 = parents.0.borrow_mut().genome.meiosis(genetic_map);

                self.genome = gamete_1.fertilize(&gamete_2);
            },
        };
        Ok(true)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Parents(Rc<RefCell< Individual>>, Rc<RefCell<Individual>>);

impl Parents{
    pub fn new(parents: (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>)) -> Parents {
        Parents(parents.0, parents.1)
    }
}

impl std::fmt::Display for Parents {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} <-> {}", 
            RefCell::borrow(&self.0).label,
            RefCell::borrow(&self.1).label
        )
    }
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
        write!(f, "tab: {: <10} label: {: <10} - parents: {: <25}", tag, self.label, parents)
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

#[derive(Debug, Clone)]
pub struct Comparison {
    label           : String,
    pair            : (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>),
    _self_comparison : bool,
}

impl Comparison {
    pub fn new(label: String, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>), self_comparison: bool) -> Comparison {
        let pair = Self::format_pair(pair);
        Comparison{label, pair, _self_comparison: self_comparison }
    }

    pub fn set_pair(&mut self, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)){
        self.pair = Self::format_pair(pair);
    }

    
    fn format_pair(pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>) {
        (Rc::clone(pair.0), Rc::clone(pair.1))
    }

    pub fn compare_genome(&self) -> (u32, u32) {
        let (mut genomewide_overlap, mut genomewide_pwd): (u32,u32)  = (0,0);
        for (chromosome1, chromosome2) in self.pair.0.borrow().genome.values().zip(self.pair.1.borrow().genome.values()){
            if ! chromosome1.is_empty() && ! chromosome2.is_empty() {
                let (chr_overlap, chr_pwd) = Self::compare_chromosome(chromosome1, chromosome2);
                genomewide_overlap += chr_overlap;
                genomewide_pwd     += chr_pwd;
            }
        }
        (genomewide_pwd, genomewide_overlap)
    }

    fn compare_chromosome(chromosome1: &Chromosome, chromosome2: &Chromosome) -> (u32, u32) {
        let (mut overlap, mut pwd): (u32,u32)  = (0,0);
        for loci in chromosome1.loci().zip(chromosome2.loci()){
            overlap+=1;
            let random_sample1 = Self::simulate_observed_reads(1, 0.0, 0.0, 0.0, loci.0.alleles());
            let random_sample2 = Self::simulate_observed_reads(1, 0.0, 0.0, 0.0, loci.1.alleles());
            if random_sample1 != random_sample2 {
                pwd+=1;
            }
        }
        (overlap, pwd)
    }

    fn simulate_observed_reads(n: u8, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64, alleles: (u8, u8)) -> Vec<u8> {
        let mut reads = Vec::new();
        let mut rng = rand::thread_rng();

        // Simulate contamination.
        for _ in 0..n {
            let chosen_base: u8 = match rng.gen::<f64>() < contam_rate {
                true  => match rng.gen::<f64>() < contam_pop_af {
                    true  => 1,  // Becomes the alternative reference allele, if contam_rate * contam_pop_af
                    false => 0,  // otherwise, pick the reference allele.
                }
                false => *[alleles.0, alleles.1].choose(&mut rng).unwrap(),
            };

            // Simulate sequencing error rate.
            let seqerror_choices: Vec<[u8; 3]> = vec![[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
            if rng.gen::<f64>() < seq_error_rate {
                let wrong_base: u8 = *seqerror_choices[chosen_base as usize].choose(&mut rng).unwrap();
                reads.push(wrong_base);
            }
            else {
                reads.push(chosen_base);
            }
        }
        reads
    }
}

//impl std::fmt::Display for Comparison {
//    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//        unsafe{
//            write!(f, "label: {: <5}\n\t- ind1: ({:?})- ind2: ({:?})", self.label, self.pair.0.borrow(), self.pair.1.borrow())
//        }
//    }
//}

#[derive(Debug, Clone)]
pub struct Pedigree {
    pub individuals: BTreeMap<String, Rc<RefCell<Individual>>>,
    comparisons: Vec<Comparison>,
    pop : Option<String>, // EUR, AFR, etc.
}

impl Pedigree {
    pub fn new() -> Pedigree {
        Pedigree { individuals: BTreeMap::new(), comparisons: Vec::new(), pop: None}
    }

    pub fn add_individual(&mut self, label: &str, parents: Option<(&String, &String)>, genome: Genome) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parents = match parents {
            None => None,
            Some((parent1, parent2)) => {
                let parent1 = self.individuals.get(parent1).ok_or(InvalidInput)?;
                let parent2 = self.individuals.get(parent2).ok_or(InvalidInput)?;
                Some((parent1, parent2))
            },
        };
        let ind = Rc::new(RefCell::new(Individual::new(label.to_owned(), parents, genome)));
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


    pub fn set_tags(&mut self, panel: &VCFPanelReader, pop: &String ) {
        self.pop = Some(pop.to_owned());
        for mut founder in self.founders_mut() {
            founder.set_tag(panel.random_sample(pop).unwrap().clone());
        }
    }

    pub fn populate_founders_vcf(&mut self, input_vcf_paths: &Vec<PathBuf>, valid_positions: &HashSet<SNPCoord>) -> Result<(), Box<dyn Error>> {
        for vcf in input_vcf_paths {
            info!("Parsing vcf file: {}", vcf.to_str().unwrap_or("None"));
            let mut vcf_reader = VCFReader::new(vcf.as_path())?;
            let pop_tag = &self.pop.as_ref().unwrap().clone();
            vcf_reader.parse_samples(self.founders_mut(), valid_positions, pop_tag)?;

        }
        Ok(())
    }
    
    pub fn reproduce(&mut self, genetic_map: &GeneticMap){
        for mut offspring in self.offsprings_mut() {
            if offspring.has_empty_genome(){
                println!("{} has empty genome.", offspring.label);
                offspring.generate_genome(genetic_map).unwrap();
            }
        }
    }

    pub fn compare_genomes(&self) {
        for comparison in &self.comparisons {
            let (pwd, overlap) = comparison.compare_genome();
            let avg_pwd = pwd as f64 / overlap as f64;
            println!("{: <20} : {: >9} - {: >9} - {: <9}", comparison.label, pwd, overlap, avg_pwd);
        }
    }
    
}

impl Default for Pedigree {
    fn default() -> Self {
        Self::new()
    }
}