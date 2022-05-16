use std::{collections::HashMap, path::{PathBuf, Path}, error::Error};

use super::PedigreeReps;
use super::pedigree::pedigree_params::ParamRateGenerator;

use crate::{
    io::{
        self,
        genotype_reader::GenotypeReader,
        vcf::{
            SampleTag,
            reader::{VCFReader, VCFPanelReader},
        }, 
    },
};

use genome::{
    SNPCoord,
    Genome,
    GeneticMap,
};

use pwd_from_stdin::{self, comparisons::Comparisons};

use rand::{Rng, prelude::ThreadRng};
use log::{info, trace, warn, debug};



pub struct Pedigrees<'a> {
    pedigrees: HashMap<String, PedigreeReps>,
    comparisons: &'a Comparisons,
    pedigree_pop: String,
    genetic_map: GeneticMap,
    previous_positions: HashMap<String, u32>,
    rng: ThreadRng,
}

/// @TODO! : Right now we're reading the pedigree definition file for each replicate. Which is highly inefficient.
///          This is because simply using clone() on a pedigree_template will merely clone() the Rc<RefCell>.
///          ==> This is NOT what we want. What we want is to create new instances of Rc<RefCell> from these values.
impl<'a> Pedigrees<'a> {



    pub fn initialize(pedigree_pop: String, comparisons: &'a Comparisons, recomb_dir: &PathBuf) -> Result<Self, Box<dyn Error>> {
        // Generate pedigree replicates for each pwd_from_stdin::Comparison.
        let pedigrees = HashMap::new();


        // --------------------- Parse input recombination maps.
        info!("Parsing genetic maps in {}", &recomb_dir.to_str().unwrap_or("None"));
        let genetic_map = GeneticMap::default().from_dir(recomb_dir)?;

        // --------------------- For each comparison, keep a record of the previously typed SNP's position.

        let mut previous_positions = HashMap::new();
        for comparison in comparisons.iter() {
            previous_positions.insert(comparison.get_pair().to_owned(), 0);
        }
        // --------------------- Initialize RNG
        let rng = rand::thread_rng();

        Ok(Pedigrees{pedigrees, comparisons, pedigree_pop, previous_positions, genetic_map, rng})
    }

    pub fn populate(&mut self, panel: &VCFPanelReader, reps: u32, genome: &Genome, pedigree_path: &Path, contam_pop: &Vec<String>, contam_num_ind: &Vec<usize>) -> Result<(), Box<dyn Error>> {
        let samples_contam_tags: Vec<Vec<SampleTag>> = panel.fetch_contaminants(contam_pop, contam_num_ind);

        for comparison in self.comparisons.iter() {
            let comparison_label = comparison.get_pair();
            let pair_indices = comparison.get_pair_indices();

            let mut pedigree_reps = PedigreeReps::with_capacity(reps as usize);
            pedigree_reps.set_contaminants(&samples_contam_tags, pair_indices);
            debug!("Contaminant set for {comparison_label}: {:#?}", pedigree_reps.contaminants);
            pedigree_reps.populate(pedigree_path, &self.pedigree_pop, panel, genome)?;
            self.pedigrees.insert(comparison_label.to_owned(), pedigree_reps);
        }
        Ok(())
    }

    pub fn set_params(&mut self, snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: &Vec<Vec<f64>>, contam_rate: &Vec<Vec<f64>>) -> Result<(), String>{

        //let samples_indices = self.comparisons.get_pairs_indices();
        for comparison in self.comparisons.iter() {
            let pair_indices = comparison.get_pair_indices();
            let pair_label = comparison.get_pair();
            let pedigree_reps = self.pedigrees.get_mut(&pair_label).ok_or(format!("Cannot access pair {pair_label}"))?;

            let mut seq_error_rate_gen = ParamRateGenerator::from_user_input(seq_error_rate, pair_indices);
            let mut contam_rate_gen = ParamRateGenerator::from_user_input(contam_rate, pair_indices);

            pedigree_reps.iter_mut()
            .for_each(|ped| {
                ped.set_params(snp_downsampling_rate, af_downsampling_rate, seq_error_rate_gen.gen_random_values(), contam_rate_gen.gen_random_values())
            });
        }
        Ok(())
    }

    fn update_pedigrees(&mut self, reader: &dyn GenotypeReader, chromosome:u8, position: u32, comparison_label: &String) -> Result<(), Box<dyn Error>>{

        // Compute the interval between current and previous position, search trough the genetic map interval tree,
        // And compute the probability of recombination.
        let previous_position    : u32 = self.previous_positions[comparison_label];
        let interval_prob_recomb : f64 = self.genetic_map.compute_recombination_prob(chromosome, previous_position, position);
        //trace!("SNP candidate for {comparison_label} - [{chromosome:<2} {position:>9}] - pop_af: {pop_af:?} - cont_af: {cont_af:?} - recomb_prop: {interval_prob_recomb:<8.6}");


        let pedigree_vec = self.pedigrees.get_mut(comparison_label).unwrap();

        // -------------------- Get contaminating population allele frequency
        //let cont_af = Some(reader.compute_local_cont_af(pedigree_vec.contaminants.as_ref())?).ok_or("Failed to retrieve contamination allele frequency")?;
        let cont_af = pedigree_vec.contaminants.as_ref().ok_or("Empty contaminants.")?.compute_local_cont_af(reader)?;

        'pedigree: for (i, pedigree) in pedigree_vec.iter_mut().enumerate() {
            // --------------------- Perform SNP downsampling if necessary
            if self.rng.gen::<f64>() < pedigree.get_params()?.snp_downsampling_rate {continue 'pedigree}

            // --------------------- Update founder alleles. Perform Allele Frequency downsampling if necessary.
            pedigree.update_founder_alleles(reader, &mut self.rng)?;

            // --------------------- Compute offspring genomes
            pedigree.compute_offspring_alleles(interval_prob_recomb, i)?;

            // --------------------- Compare genomes.
            pedigree.compare_alleles(cont_af)?;
            // --------------------- Clear genotypes before the next line!
            pedigree.clear_alleles();
        }
        Ok(())
    }

    #[allow(unused_labels)]
    pub fn pedigree_simulations_fst(&mut self, input_fst_path: &Path, maf: f64) -> Result<(), Box<dyn Error>> {
        use std::ops::Bound::*;
        // --------------------- Read and store FST index into memory.
        let mut fst_reader = io::fst::FSTReader::new(input_fst_path.to_str().unwrap());

        // Get a list of which chromosomes are found within the FST index.
        let contained_chromosomes = fst_reader.find_chromosomes();


        'comparison: for comparison in self.comparisons.iter() {
            info!("Performing simulations for : {}", comparison.get_pair());

            // Extract positions matching the FST index chromosome(s)
            // NOTE: Here we're using range(). But what we'd ideally like is drain_filter() . [Nightly only ]
            let start = SNPCoord{chromosome: *contained_chromosomes.iter().min().unwrap(), position: std::u32::MIN, reference: None, alternate: None};
            let stop =  SNPCoord{chromosome: *contained_chromosomes.iter().max().unwrap(), position: std::u32::MAX, reference: None, alternate: None};
            let relevant_positions = comparison.positions.range((Included(&start), Included(&stop)));

            // Get the number of typed positions for these chromosomes (this .clone() is quite expensive for just a fancy progress estimation...)
            let n = relevant_positions.clone().count();

            // Loop along positions.
            'coordinate: for (i, coordinate) in relevant_positions.enumerate() {


                let (chromosome, position) = (coordinate.chromosome, coordinate.position);

                // --------------------- Print progress in increments of 10%
                if (i % (n/10)) == 0 {
                    let percent = (i as f32 / (n as f32).floor()) * 100.0 ;
                    info!("{percent: >5.1}% : [{chromosome: <2} {position: >9}]");
                }

                // --------------------- Don't even begin if we know this set does not contain this chromosome
                //if ! fst_reader.contains_chr(chromosome){ continue 'coordinate }

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

                // --------------------- Parse genotype fields and start updating dynamic simulations.
                self.update_pedigrees(&fst_reader, chromosome, position, &comparison.get_pair())?;
            }
        }
        Ok(())
    }

    pub fn pedigree_simulations_vcf(&mut self, input_vcf_path: &Path, maf: f64, threads: usize) -> Result<(), Box<dyn Error>> {
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
            //let mut genotypes_filled: bool = false;
            'comparison: for comparison in self.comparisons.iter(){
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
                        let pop_af = vcf_reader.get_pop_allele_frequency(&self.pedigree_pop)?;

                        // Skip line if allele frequency is < maf
                        // [WARN]: Note that {POP}_AF entries in the INFO field relate to the REF allele frequency, not MAF.
                        //         ==> for maf = 0.05, we must also filter out alleles with AF > 0.95.
                        if pop_af < maf || pop_af > (1.0-maf) {
                                trace!("skip allele at [{chromosome: <2} {position: >9}]: pop_af: {pop_af:<8.5} --maf: {maf}");
                                vcf_reader.next_line()?;
                                continue 'line
                        }

                        // Parse genotypes. 
                        vcf_reader.fill_genotypes()?;
                    }

                    // Parse genotype fields and start updating dynamic simulations.
                    self.update_pedigrees(&vcf_reader, chromosome, position, &comparison.get_pair().to_owned())?;          
                }
            }
            // Reset line if we never parsed genotypes.
            vcf_reader.next_line()?;
        }
        Ok(())        
    }

    pub fn write_simulations(&self, output_files: &HashMap<String, String>) -> Result<(), Box<dyn Error>> {
        // --------------------- Print pedigree simulation results.
        for comparison in self.comparisons.iter() {
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

        for comparison in self.comparisons.iter() {
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