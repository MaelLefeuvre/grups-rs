use std::{collections::{BTreeMap, HashSet, HashMap}, path::{PathBuf, Path}, error::Error, ops::{Deref, DerefMut}};

use crate::io::{vcf::{
    SampleTag,
    reader::VCFReader,
    reader::VCFPanelReader
}, self, genotype_reader::GenotypeReader};

use genome::{
    SNPCoord,
    Chromosome,
    Genome,
    GeneticMap,
};

use pwd_from_stdin::{self, comparison::Comparisons};

use rand::{Rng, prelude::ThreadRng};
use rand::seq::SliceRandom;
use log::{info, trace, warn};



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

type ParentsRef<'a> = [&'a Rc<RefCell<Individual>>; 2];

impl Individual {
    pub fn new(label: String, parents: Option<ParentsRef>, genome: Genome) -> Individual {
        let parents = parents.map(Self::format_parents);
        Individual {tag: None, label, parents, genome, strands: None, currently_recombining: [false, false], alleles: None}
    }

    pub fn set_alleles(&mut self, alleles: [u8; 2]) {
        self.alleles = Some(alleles);
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
        self.parents.as_ref().map(|parents| (RefCell::borrow(&parents[0]), RefCell::borrow(&parents[1])))
    }

    pub fn get_parents_mut(&self) -> Option<[RefMut<Individual>; 2]> {
        self.parents.as_ref().map(|parents| [RefCell::borrow_mut(&parents[0]), RefCell::borrow_mut(&parents[1])])
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

    pub fn has_empty_genome(&self) -> bool {
        self.genome.is_empty()
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


    pub fn generate_genome(&mut self, genetic_map: &GeneticMap) -> Result<bool, Box<dyn Error>>{
        if ! self.has_empty_genome() {
            trace!("{} Genome already generated.", self.label);
            return Ok(false)
        }

        trace!("Generating {} genome.", self.label);
        match &self.parents {
            None => panic!("Cannot generate genome, as parents are missing."),
            Some(parents) => {
                for parent in parents.iter() {
                    if parent.borrow().has_empty_genome() {
                        trace!(" -> Parent {} has empty genome. Generating...", parent.borrow().label);
                        parent.borrow_mut().generate_genome(genetic_map).unwrap();
                    }
                    else {
                        trace!("-> Parent {} has genome.", parent.borrow().label);
                    }
                }
                let gamete_1 = parents[0].borrow_mut().genome.meiosis(genetic_map);
                let gamete_2 = parents[1].borrow_mut().genome.meiosis(genetic_map);

                self.genome = gamete_1.fertilize(&gamete_2);
            },
        };
        Ok(true)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Parents([Rc<RefCell<Individual>>; 2]);

impl Parents{
    pub fn new(parents: [Rc<RefCell<Individual>>; 2]) -> Parents {
        Parents(parents)
    }
}

impl Deref for Parents {
    type Target = [Rc<RefCell<Individual>>; 2];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::fmt::Display for Parents {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} <-> {}", 
            RefCell::borrow(&self[0]).label,
            RefCell::borrow(&self[1]).label
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

#[derive(Debug, Clone)]
pub struct PedComparison {
    label            : String,
    pair             : (Rc<RefCell<Individual>>, Rc<RefCell<Individual>>),
    pwd              : u32,
    overlap          : u32,
    _self_comparison : bool,
}

impl PedComparison {
    pub fn new(label: String, pair: (&Rc<RefCell<Individual>>, &Rc<RefCell<Individual>>), self_comparison: bool) -> PedComparison {
        let pair = Self::format_pair(pair);
        PedComparison{label, pair, _self_comparison: self_comparison, pwd: 0, overlap: 0 }
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

impl std::fmt::Display for PedComparison {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let default_tag=SampleTag::new("None", None);
        write!(f, " {: <20} {: <8} {: <8} {: <8} {: <8} {: >9} {: >9} {: <12.6}",
            self.label,
            self.pair.0.borrow().label,
            self.pair.1.borrow().label,
            self.pair.0.borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pair.1.borrow().tag.as_ref().unwrap_or(&default_tag).id(),
            self.pwd,
            self.overlap,
            self.get_avg_pwd()
        )
    }
}

#[derive(Debug, Clone)]
struct PedComparisons(Vec<PedComparison>);

impl Deref for PedComparisons {
    type Target = Vec<PedComparison>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for PedComparisons {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


impl std::fmt::Display for PedComparisons {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0.iter().fold(Ok(()), |result, comparison| {
            result.and_then(|_| writeln!(f, "{}", comparison))
        })
    }
}

impl PedComparisons {
    pub fn new() -> Self {
        Self(Vec::new())
    }
}

#[derive(Debug, Clone)]
pub struct Pedigree {
    pub individuals: BTreeMap<String, Rc<RefCell<Individual>>>,
    comparisons: PedComparisons,
    pop : Option<String>, // EUR, AFR, etc.
}

impl Pedigree {
    pub fn new() -> Pedigree {
        Pedigree { individuals: BTreeMap::new(), comparisons: PedComparisons::new(), pop: None}
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
                Some([parent1, parent2])
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
        self.comparisons.push(PedComparison::new(label.to_owned(), (pair0, pair1), pair0==pair1));
        Ok(())
    }

    pub fn set_relationship(&mut self, ind: &String, parents: (&String, &String)) -> std::io::Result<()>{
        use std::io::ErrorKind::InvalidInput;
        let parent0 = &self.individuals.get(parents.0).ok_or(InvalidInput)?.clone();
        let parent1 = &self.individuals.get(parents.1).ok_or(InvalidInput)?.clone();

        self.individuals.get_mut(ind)
            .ok_or(InvalidInput)?
            .borrow_mut()
            .set_parents([parent0, parent1]);
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
            let mut vcf_reader = VCFReader::new(vcf.as_path(), 0)?;
            let pop_tag = &self.pop.as_ref().unwrap().clone();
            vcf_reader.parse_samples(self.founders_mut(), valid_positions, pop_tag)?;

        }
        Ok(())
    }
    
    pub fn reproduce(&mut self, genetic_map: &GeneticMap){
        for mut offspring in self.offsprings_mut() {
            if offspring.has_empty_genome(){
                trace!("{} has empty genome.", offspring.label);
                offspring.generate_genome(genetic_map).unwrap();
            }
        }
    }

    pub fn compare_genomes(&self) {
        for comparison in self.comparisons.iter() {
            let (pwd, overlap) = comparison.compare_genome();
            let avg_pwd = pwd as f64 / overlap as f64;
            println!("{: <20} : {: >9} - {: >9} - {: <9}", comparison.label, pwd, overlap, avg_pwd);
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

pub struct PedigreeReps(Vec<Pedigree>);

impl PedigreeReps{
    pub fn with_capacity(n: usize) -> Self {
        Self(Vec::with_capacity(n))
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
        &self.0
    }
}

impl DerefMut for PedigreeReps {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl std::fmt::Display for PedigreeReps {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0.iter().enumerate().fold(Ok(()), |_, (idx, pedigree)| {
            pedigree.comparisons.iter().fold(Ok(()), |result, comparison| {
                result.and_then(|_| writeln!(f, "{idx} {}", comparison))
            })
        })
    }
}


pub struct Pedigrees<'a> {
    pedigrees: HashMap<String, PedigreeReps>,
    comparisons: &'a Comparisons,
    pedigree_pop: String,
    genetic_map: GeneticMap,
    rng: ThreadRng,
}

impl<'a> Pedigrees<'a> {
    pub fn new(pedigree_path: &Path, reps: u32, pedigree_pop: String, comparisons: &'a Comparisons, panel: &VCFPanelReader, genome: &Genome, recomb_dir: &PathBuf) -> Result<Self, Box<dyn Error>> {

        // Create template pedigree from definition file.
        let template_pedigree= io::pedigree::pedigree_parser(pedigree_path, genome).unwrap();

        // Generate pedigree replicates for each pwd_from_stdin::Comparison.
        let mut pedigrees = HashMap::new();
        for comparison in comparisons.get() {
            let comparison_label = comparison.get_pair();
            pedigrees.insert(comparison_label.to_owned(), PedigreeReps::with_capacity(reps as usize));
            for _ in 0..reps {
                let mut new_pedigree = template_pedigree.clone();
                new_pedigree.set_tags(panel, &pedigree_pop);
                new_pedigree.assign_offspring_strands()?;
                pedigrees.get_mut(&comparison_label).unwrap().push(new_pedigree);
            }
        }

        // --------------------- Parse input recombination maps.
        info!("Parsing genetic maps in {}", &recomb_dir.to_str().unwrap_or("None"));
        let genetic_map = GeneticMap::default().from_dir(recomb_dir)?;

        // --------------------- Initialize RNG
        let rng = rand::thread_rng();

        Ok(Pedigrees{pedigrees, comparisons, pedigree_pop, genetic_map, rng})
    }


    fn update_pedigrees(&mut self, reader: &dyn GenotypeReader, comparison_label: &String, snp_downsampling_rate: f64, af_downsampling_rate: f64, interval_prob_recomb: f64, contam_rate: f64, cont_af: f64, seq_error_rate: f64) -> Result<(), Box<dyn Error>>{
        let pedigree_vec = self.pedigrees.get_mut(comparison_label).unwrap();
        'pedigree: for (i, pedigree) in pedigree_vec.iter_mut().enumerate() {
            // --------------------- Perform SNP downsampling if necessary
            if self.rng.gen::<f64>() < snp_downsampling_rate {continue 'pedigree}

            // --------------------- Update founder alleles. Perform Allele Frequency downsampling if necessary.
            for mut founder in pedigree.founders_mut() {
                founder.alleles = match self.rng.gen::<f64>() < af_downsampling_rate { 
                    false => reader.get_alleles(founder.get_tag().unwrap()),
                    true  => Some([0, 0]),
                };
            }

            // --------------------- Compute offspring genomes
            for mut offspring in pedigree.offsprings_mut() {
                offspring.assign_alleles(interval_prob_recomb, i)?;
            }

            // --------------------- Compare genomes.
            for comparison in &mut pedigree.comparisons.iter_mut() {
                comparison.compare_alleles(contam_rate, cont_af, seq_error_rate);
            }
            // --------------------- Clear genotypes before the next line!
            pedigree.clear_alleles();
        }
        Ok(())
    }

    #[allow(unused_labels)]
    pub fn pedigree_simulations_fst(&mut self, input_fst_path: &Path, contam_rate: f64, contam_ind_ids: &[&SampleTag], seq_error_rate: f64, af_downsampling_rate: f64, snp_downsampling_rate: f64, maf: f64) -> Result<(), Box<dyn Error>> {
        // --------------------- For each comparison, keep a record of the previously typed SNP's position.       // DUPLICATE
        let mut previous_positions = HashMap::new();                                        // DUPLICATE
        for comparison in self.comparisons.get() {                                                   // DUPLICATE
            previous_positions.insert(comparison.get_pair().to_owned(), 0);                                  // DUPLICATE
        }

        // --------------------- Read and store FST index into memory.
        let mut fst_reader = io::fst::FSTReader::new(input_fst_path.to_str().unwrap());

        'comparison: for comparison in self.comparisons.get() {
            info!("Performing simulations for : {}", comparison.get_pair());
            'coordinate: for (i, coordinate) in comparison.positions.iter().enumerate() {

                let (chromosome, position) = (coordinate.chromosome, coordinate.position);

                // --------------------- Print progress in increments of 10%
                if (i % (comparison.positions.len()/10)) == 0 {
                    let percent = (i as f32 / (comparison.positions.len() as f32).floor()) * 100.0 ;
                    info!("{percent: >5.1}% : [{chromosome: <2} {position: >9}]");
                }

                // --------------------- Don't even begin if we know this set does not contain this chromosome
                if ! fst_reader.contains_chr(chromosome){
                    continue 'coordinate
                }

                // --------------------- Search through the FST index for the genotypes and pop frequencies at this coordinate.
                fst_reader.clear_buffers();
                fst_reader.search_coordinate_genotypes(chromosome, position);
                fst_reader.search_coordinate_frequencies(chromosome, position);

                // --------------------- If genotypes is empty, then this position is missing within the index... Skip ahead.
                if ! fst_reader.has_genotypes() {
                    warn!("Missing coordinate in fst index: [{chromosome: <2} {position: >9}]");
                    continue 'coordinate
                }

                // --------------------- Skip line if allele frequency is < maf
                let pop_af = fst_reader.get_pop_allele_frequency(&self.pedigree_pop)?;
                if pop_af < maf || pop_af > (1.0-maf) { 
                    trace!("skip allele at [{chromosome: <2} {position: >9}]: pop_af: {pop_af:<8.5} --maf: {maf}");
                    continue 'coordinate
                }

                // -------------------- Get contaminating population allele frequency
                let cont_af = Some(fst_reader.compute_local_cont_af(contam_ind_ids)?);


                // --------------------- Compute the probability of recombination using our genetic map.
                let comparison_label = comparison.get_pair().to_owned();
                let previous_position = previous_positions[&comparison_label];
                let interval_prob_recomb: f64 = self.genetic_map.compute_recombination_prob(chromosome, previous_position, position);
                trace!("SNP candidate for {comparison_label} - [{chromosome:<2} {position:>9}] - pop_af: {pop_af:?} - cont_af: {cont_af:?} - recomb_prop: {interval_prob_recomb:<8.6}");
                // --------------------- Parse genotype fields and start updating dynamic simulations.
                self.update_pedigrees(&fst_reader, &comparison_label, snp_downsampling_rate, af_downsampling_rate, interval_prob_recomb, contam_rate, cont_af.unwrap(), seq_error_rate)?;
            }
        }
        Ok(())
    }

    pub fn pedigree_simulations_vcf(&mut self, input_vcf_path: &Path, contam_rate: f64, contam_ind_ids: &[&SampleTag], seq_error_rate: f64, af_downsampling_rate: f64, snp_downsampling_rate: f64, maf: f64, threads: usize) -> Result<(), Box<dyn Error>> {
        // --------------------- For each comparison, keep a record of the previously typed SNP's position.       // DUPLICATE
        let mut previous_positions = HashMap::new();                                        // DUPLICATE
        for comparison in self.comparisons.get() {                                                   // DUPLICATE
            previous_positions.insert(comparison.get_pair().to_owned(), 0);                                  // DUPLICATE
        }

        // --------------------- Read VCF File line by line.
        let mut i = 0;
        let mut vcf_reader = VCFReader::new(input_vcf_path, threads)?;
        'line: while vcf_reader.has_data_left()? {

            // Get current chromosome and position.
            let (chromosome, position) = vcf_reader.parse_coordinate()?; // 1

            // --------------------- Print progress in increments of 50 000 lines.
            if i % 50_000 == 0 { info!(" {i: >9}: [{chromosome: <2} {position: >9}]");}
            i+=1;

            // --------------------- Loop across comparisons. 
            let mut pop_af : Option<f64> = None;
            let mut cont_af: Option<f64> = None;
            //let mut genotypes_filled: bool = false;
            'comparison: for comparison in self.comparisons.get(){
                // --------------------- Skip if the current position is not a valid candidate.
                if ! comparison.positions.contains(&SNPCoord{chromosome, position, reference: None, alternate: None}){
                    continue 'comparison
                } else {
                    if ! vcf_reader.genotypes_filled {
                        // Go to INFO field.
                        vcf_reader.parse_info_field()?;

                        // Check if this coordinate is a biallelic SNP and skip line if not.
                        if vcf_reader.is_multiallelic() || ! vcf_reader.is_snp() {
                            vcf_reader.next_line()?;
                            continue 'line
                        }

                        // Extract population allele frequency.
                        pop_af = Some(vcf_reader.get_pop_allele_frequency(&self.pedigree_pop)?);

                        // Skip line if allele frequency is < maf
                        // [WARN]: Note that {POP}_AF entries in the INFO field relate to the REF allele frequency, not MAF.
                        //         ==> for maf = 0.05, we must also filter out alleles with AF > 0.95.
                        match pop_af {
                            Some(af) => if af < maf || af > (1.0-maf) {
                                trace!("skip allele at [{chromosome: <2} {position: >9}]: pop_af: {af:<8.5} --maf: {maf}");
                                vcf_reader.next_line()?;
                                continue 'line
                            },
                            None => panic!("Empty population allele frequency!")
                        }

                        // Parse genotypes. 
                        vcf_reader.fill_genotypes()?;
                        //genotypes_filled = true;
    
                        // Compute contaminating pop allele frequency
                        cont_af = Some(vcf_reader.compute_local_cont_af(contam_ind_ids)?);
                    }

                    // Compute the interval between current and previous position, search trough the genetic map interval tree,
                    // And compute the probability of recombination.
                    let comparison_label = comparison.get_pair().to_owned();
                    let previous_position = previous_positions[&comparison_label];
                    let interval_prob_recomb: f64 = self.genetic_map.compute_recombination_prob(chromosome, previous_position, position);
                    trace!("SNP candidate for {comparison_label} - [{chromosome:<2} {position:>9}] - pop_af: {pop_af:?} - cont_af: {cont_af:?} - recomb_prop: {interval_prob_recomb:<8.6}");

                    // Parse genotype fields and start updating dynamic simulations.
                    self.update_pedigrees(&vcf_reader, &comparison_label, snp_downsampling_rate, af_downsampling_rate, interval_prob_recomb, contam_rate, cont_af.unwrap(), seq_error_rate)?;          
                }
            }
            // Reset line if we never parsed genotypes.
            vcf_reader.next_line()?;
        }
        Ok(())        
    }

    pub fn write_simulations(&self, output_files: &HashMap<String, String>) -> Result<(), Box<dyn Error>> {
        // --------------------- Print pedigree simulation results.
        for comparison in self.comparisons.get() {
            let comparison_label = comparison.get_pair();
            let pedigree_vec = self.pedigrees.get(&comparison_label).unwrap();

            let mut writer = pwd_from_stdin::io::Writer::new(Some(output_files[&comparison_label].clone()))?;
            writer.write_iter(vec![&pedigree_vec])?;

            info!("\n--------------------- {comparison_label}\n{}", pedigree_vec);
        }
        Ok(())
    }

    pub fn compute_results(&self, output_file: &String) -> Result<(), Box<dyn Error>> {
        let mut simulations_results = Vec::with_capacity(self.comparisons.len());

        for comparison in self.comparisons.get() {
            let comparison_label    = comparison.get_pair();
            let pedigree_vec = self.pedigrees.get(&comparison_label).unwrap();
    
            let sum_simulated_pwds = pedigree_vec.compute_sum_simulated_pwds();
    

            let observed_avg_pwd        : f64 = comparison.get_avg_pwd();
            let mut min_z_score         : f64 = f64::MAX;
            let mut most_likely_avg_pwd : f64 = 0.0;
            let mut most_likely_rel     : String = "None".to_string();

            for (scenario, simulated_sum_avg_pwd) in sum_simulated_pwds.iter() {
                let avg_avg_pwd: f64 = simulated_sum_avg_pwd/pedigree_vec.len() as f64;
                let scenario_z_score = (avg_avg_pwd - observed_avg_pwd).abs();
                if  scenario_z_score < min_z_score {
                    min_z_score =  scenario_z_score;
                    most_likely_rel = scenario.to_owned();
                    most_likely_avg_pwd = avg_avg_pwd;
                }
            }
    
            let min_z_score = if min_z_score == f64::MAX {f64::NAN} else {min_z_score};
            let pair_name = comparison.get_pair();

            let simulation_result = format!("{pair_name: <20} {most_likely_rel: <20} {observed_avg_pwd: >8.6} {most_likely_avg_pwd: >8.6} {min_z_score: <8.6}");
            println!("{simulation_result}");
            simulations_results.push(simulation_result);
        }


        let mut writer = pwd_from_stdin::io::Writer::new(Some(output_file.clone()))?;
        writer.write_iter(simulations_results)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_pedigree() -> std::io::Result<Pedigree> {
        let mut pedigree = Pedigree::new();
        pedigree.add_individual("father", None, Genome::default())?;
        pedigree.add_individual("mother", None, Genome::default())?;
        pedigree.add_individual("offspr", Some((&"father".to_string(), &"mother".to_string())), Genome::default())?;

        let mut father = pedigree.get_mutind(&"father".to_string()).expect("Cannot extract father");
        father.set_alleles([0, 1]);
        drop(father);

        let mut mother = pedigree.get_mutind(&"mother".to_string()).expect("Cannot extract mother");
        mother.set_alleles([1, 0]);
        drop(mother);


        Ok(pedigree)
    }

    #[test]
    #[should_panic]
    fn meiosis_assign_alleles_empty_strands(){
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
    }

    #[test]
    fn meiosis_assign_alleles_filled_strands(){
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,0]);

        let output = offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(output, true);

        let output = offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(output, false);
    }

    #[test]
    fn meiosis_check_strands_00() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");

        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,0]);
        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([0, 1]))
    }

    #[test]
    fn meiosis_check_strands_01() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([1,1]);

        offspr.assign_alleles(0.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 0]));
        assert_eq!(offspr.currently_recombining, [false, false]);

    }

    #[test]
    fn meiosis_check_recombination() {
        let mut pedigree = test_pedigree().expect("Cannot generate test pedigree");
        let mut offspr = pedigree.get_mutind(&"offspr".to_string()).expect("Cannot extract offspr");
        offspr.strands= Some([0,1]);
        offspr.assign_alleles(1.0, 0).expect("Failed to assign alleles");
        assert_eq!(offspr.alleles, Some([1, 1]));
        assert_eq!(offspr.currently_recombining, [true, true]);

    }



}