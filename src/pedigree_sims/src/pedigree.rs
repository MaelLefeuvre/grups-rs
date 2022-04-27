use std::{collections::{BTreeMap, HashSet, HashMap}, path::PathBuf, error::Error};
use crate::io::{SampleTag, VCFReader, VCFAsyncReader, VCFPanelReader};
use log::{info, debug};
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
    strands    : Option<[usize; 2]>,
    currently_recombining: [bool; 2],
    pub alleles: Option<[u8; 2]>,
}

type ParentsRef<'a> = (&'a Rc<RefCell<Individual>>, &'a Rc<RefCell<Individual>>);

impl Individual {
    pub fn new(label: String, parents: Option<ParentsRef>, genome: Genome) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label, parents, genome, strands: None, currently_recombining: [false, false], alleles: None}
    }

    pub fn alleles_tuple(&self) -> (u8, u8) {
        match self.alleles {
            None => panic!("Missing alleles."),
            Some(alleles) => (alleles[0], alleles[1])
        }
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

    pub fn clear_alleles(&mut self){
        self.alleles = None
    }
    pub fn assign_alleles(&mut self, recombination_prob: f64) -> Result<bool, Box<dyn Error>> {
        if self.alleles != None {
            //println!("{} Genome already generated.", self.label);
            return Ok(false)
        }

        match &self.parents {
            None => panic!("Cannot generate genome, as parents are missing."),
            Some(parents) => {
                let mut i = 0;
                let mut rng = rand::thread_rng();
                for parent in [&parents.0, &parents.1] {

                    //Assign parent genome if not previously generated.
                    if parent.borrow().alleles == None {
                        parent.borrow_mut().assign_alleles(recombination_prob).unwrap();
                    }

                    // Check if recombination occured for each parent and update counters if so.

                    if rng.gen::<f64>() < recombination_prob {
                        //println!("Switch!");
                        self.currently_recombining[i] = ! self.currently_recombining[i];
                    }
                    i+=1;
                }

                // Assign alleles.
                let haplo_0 = parents.0.borrow_mut().meiosis(self.strands.unwrap()[0], self.currently_recombining[0]);
                let haplo_1 = parents.1.borrow_mut().meiosis(self.strands.unwrap()[1], self.currently_recombining[1]);
                self.alleles = Some([haplo_0, haplo_1]);
            }
        }
        Ok(true)
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
                let gamete_2 = parents.1.borrow_mut().genome.meiosis(genetic_map);

                self.genome = gamete_1.fertilize(&gamete_2);
            },
        };
        Ok(true)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Parents(Rc<RefCell<Individual>>, Rc<RefCell<Individual>>);

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
    label            : String,
    pair             : (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>),
    pwd              : u32,
    overlap          : u32,
    _self_comparison : bool,
}

impl Comparison {
    pub fn new(label: String, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>), self_comparison: bool) -> Comparison {
        let pair = Self::format_pair(pair);
        Comparison{label, pair, _self_comparison: self_comparison, pwd: 0, overlap: 0 }
    }

    pub fn add_pwd(&mut self){
        self.pwd += 1;
    }

    pub fn add_overlap(&mut self){
        self.overlap +=1;
    }

    pub fn get_avg_pwd(&self) -> f64 {
        self.pwd as f64 / self.overlap as f64
    }

    pub fn set_pair(&mut self, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)){
        self.pair = Self::format_pair(pair);
    }

    
    fn format_pair(pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>)) -> (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>) {
        (Rc::clone(pair.0), Rc::clone(pair.1))
    }

    pub fn compare_alleles(&mut self, contam_rate: f64, contam_pop_af: f64, seq_error_rate: f64) {
        self.add_overlap();
        let random_sample0 = Self::simulate_observed_reads(1, contam_rate, contam_pop_af, seq_error_rate, self.pair.0.borrow().alleles_tuple());
        let random_sample1 = Self::simulate_observed_reads(1, contam_rate, contam_pop_af, seq_error_rate, self.pair.1.borrow().alleles_tuple());
        if random_sample0 != random_sample1 {
            self.add_pwd();
        }
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

    pub fn assign_offspring_strands(&mut self) -> Result<(), String> {
        for mut offspring in self.offsprings_mut() {
            offspring.assign_strands()?;
        }
        Ok(())
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
    
    pub fn print_results(&self, idx: usize) {
        for comparison in &self.comparisons {
            let label = comparison.label.to_owned();
            let pwd = comparison.pwd;
            let overlap = comparison.overlap;
            let avg_pwd = pwd as f64 / overlap as f64;
            println!("{idx} {label: <20}: {pwd: >9} {overlap: >9} {avg_pwd: <12.6}")
        }
    }

    pub fn clear_alleles(&mut self){
        for ind in self.individuals.values_mut(){
            ind.borrow_mut().clear_alleles()
        }
    }

}

impl Default for Pedigree {
    fn default() -> Self {
        Self::new()
    }
}

pub async fn pedigree_simulations(pedigrees: &mut Vec<Pedigree>, input_vcf_path: &PathBuf, valid_positions: &HashSet<SNPCoord>, pop: &String, genetic_map: &GeneticMap, threads: usize) -> Result<(), Box<dyn Error>>{
    let mut vcf_reader = VCFAsyncReader::new(input_vcf_path.as_path(), threads).await?;
    let mut previous_position = 0;
    let mut i = 0;
    while vcf_reader.has_data_left().await? {

        // Get current chromosome and position.
        let chromosome : u8  = vcf_reader.next_field().await?.parse()?; // 1
        let position   : u32 = vcf_reader.next_field().await?.parse()?; // 2
        
        if i % 50_000 == 0 {
            debug!("{i: >9} {chromosome: >2} {position: >9}");
        }
        i+=1;

        // Check if the current position is a valid candidate. Skip if not.
        if ! valid_positions.contains(&SNPCoord{chromosome, position, reference: None, alternate: None}){
            vcf_reader.skip_line().await?;
            continue
        }

        // Go to INFO field.
        vcf_reader.skip(5).await?;                                      // 6 
        let info = vcf_reader.next_field().await?.split(';').collect::<Vec<&str>>();

        // Check if Bi-Allelic and skip line if not.
        if info.iter().any(|&field| field == "MULTI_ALLELIC") {
            vcf_reader.skip_line().await?;
            continue
        }

        // Get the mutation type from INFO...
        let vtype = info.iter()
            .find(|&&field| field.starts_with("VT=")).unwrap()
            .split('=')
            .collect::<Vec<&str>>()[1];

        // ...and skip if this is not an snp.
        if vtype != "SNP" {
            vcf_reader.skip_line().await?;
            continue
        }

        // Extract population allele frequency.
        let pop_af = info.iter()
            .find(|&&field| field.starts_with(&format!("{pop}_AF")))
            .unwrap()
            .split('=')
            .collect::<Vec<&str>>()[1]
            .parse::<f64>()
            .unwrap();

        // Compute the interval between current and previous position, search trough the genetic map interval tree,
        // And compute the probability of recombination.
        let mut interval_prob_recomb: f64 = 0.0;
        for recombination_range in genetic_map[&chromosome].find(previous_position, position) {
            let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
            let real_stop  = if position          > recombination_range.stop  {recombination_range.stop } else {position         };

            interval_prob_recomb += recombination_range.val.prob() * (real_stop as f64 - real_start as f64 + 1.0);
        }

        //println!("{interval_prob_recomb}");

        // Parse genotype fields and start updating dynamic simulations.
        vcf_reader.fill_genotypes().await?;
        for pedigree in pedigrees.iter_mut() {

            // Update founder alleles.
            for mut founder in pedigree.founders_mut() {
                founder.alleles = vcf_reader.get_alleles2(founder.get_tag().unwrap().idx())?;
            }

            //Compute offspring genomes
            for mut offspring in pedigree.offsprings_mut() {
                offspring.assign_alleles(interval_prob_recomb)?;
            }

            // Compare genomes.
            for comparison in &mut pedigree.comparisons {
                comparison.compare_alleles(0.0, 0.0, 0.0);
            }

            // Clear genotypes before the next line!
            pedigree.clear_alleles();
        }
    previous_position = position;

    }
    Ok(())
}

pub fn compute_results(pedigrees: &Vec<Pedigree>, comparison: &pwd_from_stdin::comparison::Comparison) {
    let mut avg_simulated_pwd = HashMap::new();
    for pedigree in pedigrees.iter() {
        for comparison in pedigree.comparisons.iter() {
            //intervals.entry(chr).or_insert_with(Vec::new).push(interval);
            *avg_simulated_pwd.entry(comparison.label.to_owned()).or_insert(0.0) += comparison.get_avg_pwd()
        }
    }

    let mut min_z_score = f64::MAX;
    let mut most_likely_rel = "None".to_string();
    let mut most_likely_avg_pwd = 0.0;
    let observed_avg_pwd = comparison.get_avg_pwd();
    for (scenario, simulated_sum_avg_pwd) in avg_simulated_pwd.iter() {
        let avg_avg_pwd = simulated_sum_avg_pwd/pedigrees.len() as f64;
        let scenario_z_score = (avg_avg_pwd - observed_avg_pwd).abs();
        if  scenario_z_score < min_z_score {
            min_z_score =  scenario_z_score;
            most_likely_rel = scenario.to_owned();
            most_likely_avg_pwd = avg_avg_pwd;
        }
    }

    println!("{: <20} {: <20} {: >8.6} {: >8.6} {: <8.6}", comparison.get_pair(), most_likely_rel, observed_avg_pwd, most_likely_avg_pwd, min_z_score);
}