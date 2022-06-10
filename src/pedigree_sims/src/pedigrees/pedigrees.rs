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

use genome::GeneticMap;

use pwd_from_stdin::{self, comparisons::Comparisons};

use rand::{Rng, prelude::ThreadRng};
use log::{info, trace, warn, debug};

use pwd_from_stdin::comparisons::pwd::Coordinate;

pub struct Pedigrees {
    pedigrees: HashMap<String, PedigreeReps>,
    pedigree_pop: String,
    genetic_map: GeneticMap,
    previous_positions: HashMap<String, u32>,
    rng: ThreadRng,
}

/// @TODO! : Right now we're reading the pedigree definition file for each replicate. Which is highly inefficient.
///          This is because simply using clone() on a pedigree_template will merely clone() the Rc<RefCell>.
///          ==> This is NOT what we want. What we want is to create new instances of Rc<RefCell> from these values.
impl Pedigrees {
    pub fn initialize(pedigree_pop: String, comparisons: &Comparisons, recomb_dir: &PathBuf) -> Result<Self, Box<dyn Error>> {
        // Generate pedigree replicates for each pwd_from_stdin::Comparison.
        let pedigrees = HashMap::new();


        // --------------------- Parse input recombination maps.
        info!("Parsing genetic maps in {}", &recomb_dir.to_str().unwrap_or("None"));
        let genetic_map = GeneticMap::default().from_dir(recomb_dir)?;

        // --------------------- For each comparison, keep a record of the previously typed SNP's position.

        let mut previous_positions = HashMap::new();
        for comparison in comparisons.iter() {
            let key = comparison.get_pair();
            previous_positions.insert(key.to_owned(), 0);
        }
        // --------------------- Initialize RNG
        let rng = rand::thread_rng();

        Ok(Pedigrees{pedigrees, pedigree_pop, previous_positions, genetic_map, rng})
    }

    pub fn populate(&mut self, comparisons: &Comparisons,  panel: &VCFPanelReader, reps: u32, pedigree_path: &Path, contam_pop: &Vec<String>, contam_num_ind: &Vec<usize>) -> Result<(), Box<dyn Error>> {
        let samples_contam_tags: Vec<Vec<SampleTag>> = panel.fetch_contaminants(contam_pop, contam_num_ind);

        for comparison in comparisons.iter() {
            let comparison_label = comparison.get_pair();
            let pair_indices = comparison.get_pair_indices();

            let mut pedigree_reps = PedigreeReps::with_capacity(reps as usize);
            pedigree_reps.set_contaminants(&samples_contam_tags, pair_indices);
            debug!("Contaminant set for {comparison_label}: {:#?}", pedigree_reps.contaminants);
            pedigree_reps.populate(pedigree_path, &self.pedigree_pop, panel)?;
            self.pedigrees.insert(comparison_label.to_owned(), pedigree_reps);
        }
        Ok(())
    }

    pub fn set_params(&mut self, comparisons: &Comparisons, snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: &Option<Vec<Vec<f64>>>, contam_rate: &Vec<Vec<f64>>) -> Result<(), String>{

        //let samples_indices = self.comparisons.get_pairs_indices();
        for comparison in comparisons.iter() {
            let pair_indices = comparison.get_pair_indices();
            let pair_label = comparison.get_pair();
            let pedigree_reps = self.pedigrees.get_mut(&pair_label).ok_or(format!("Cannot access pair {pair_label}"))?;

            let mut seq_error_rate_gen = match seq_error_rate {
                Some(seq_errors_vec) => Some(ParamRateGenerator::from_user_input(seq_errors_vec, pair_indices)),
                None => None,
            };
            let mut contam_rate_gen = ParamRateGenerator::from_user_input(contam_rate, pair_indices);

            pedigree_reps.iter_mut()
            .for_each(|ped| {
                let seq_error_rate = match &mut seq_error_rate_gen {
                    Some(generator) => Some(generator.gen_random_values()),
                    None => None
                };

                ped.set_params(snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate_gen.gen_random_values())
            });
        }
        Ok(())
    }

    fn update_pedigrees(&mut self, reader: &dyn GenotypeReader, chromosome:u8, position: u32, comparison_label: &String, pileup_error_probs: &[f64;2]) -> Result<(), Box<dyn Error>>{

        // Compute the interval between current and previous position, search trough the genetic map interval tree,
        // And compute the probability of recombination.
        let previous_position    : u32 = self.previous_positions[comparison_label];
        let interval_prob_recomb : f64 = self.genetic_map.compute_recombination_prob(chromosome, previous_position, position);
        //trace!("SNP candidate for {comparison_label} - [{chromosome:<2} {position:>9}] - pop_af: {pop_af:?} - cont_af: {cont_af:?} - recomb_prop: {interval_prob_recomb:<8.6}");


        let pedigree_vec = self.pedigrees.get_mut(comparison_label).unwrap();

        // -------------------- Get contaminating population allele frequency
        let cont_af = pedigree_vec.contaminants.as_ref().ok_or("Empty contaminants.")?.compute_local_cont_af(reader)?;

        'pedigree: for (i, pedigree) in pedigree_vec.iter_mut().enumerate() {
            // --------------------- Perform SNP downsampling if necessary
            if self.rng.gen::<f64>() < pedigree.get_params()?.snp_downsampling_rate {continue 'pedigree}

            // --------------------- Update founder alleles. Perform Allele Frequency downsampling if necessary.
            pedigree.update_founder_alleles(reader, &mut self.rng)?;

            // --------------------- Compute offspring genomes
            pedigree.compute_offspring_alleles(interval_prob_recomb, i)?;

            // --------------------- Compare genomes.
            pedigree.compare_alleles(cont_af, &pileup_error_probs)?;
            // --------------------- Clear genotypes before the next line!
            pedigree.clear_alleles();
        }
        Ok(())
    }

    #[allow(unused_labels)]
    pub fn pedigree_simulations_fst(&mut self, comparisons: &mut Comparisons, input_fst_path: &Path, maf: f64) -> Result<(), Box<dyn Error>> {
        // --------------------- Read and store FST index into memory.
        let mut fst_reader = io::fst::FSTReader::new(input_fst_path.to_str().unwrap());

        //// Get a list of which chromosomes are found within the FST index.
        let contained_chromosomes = fst_reader.find_chromosomes();

        // Keep a record of positions that should get filtered out after maf < treshold.
        let mut positions_to_delete = HashMap::new();
        'comparison: for comparison in comparisons.iter() {
            info!("Performing simulations for : {}", comparison.get_pair());

            //// Extract positions matching the FST index chromosome(s)
            //// NOTE: Here we're using range(). But what we'd ideally like is drain_filter() . [Nightly only]
            let start = Coordinate::new(*contained_chromosomes.iter().min().unwrap(), std::u32::MIN);
            let stop  = Coordinate::new(*contained_chromosomes.iter().max().unwrap(), std::u32::MAX);
            let relevant_positions = comparison.positions.range(start..stop);
            //
            // Get the number of typed positions for these chromosomes (this .clone() is quite expensive for just a fancy progress estimation...)
            let n = relevant_positions.clone().count();

            // Loop along positions.
            let key = comparison.get_pair();
            'coordinate: for (i, pairwise_diff) in relevant_positions.enumerate() {                
                let (chromosome, position) = (pairwise_diff.coordinate.chromosome, pairwise_diff.coordinate.position);

                //// --------------------- Print progress in increments of 10%
                if (i % ((n/10)+1)) == 0 {
                    let percent = (i as f32 / (n as f32).floor()) * 100.0 ;
                    info!("{percent: >5.1}% : [{chromosome: <2} {position: >9}]");
                }

                // --------------------- Search through the FST index for the genotypes and pop frequencies at this coordinate.
                fst_reader.clear_buffers();
                fst_reader.search_coordinate_genotypes(chromosome, position);
                fst_reader.search_coordinate_frequencies(chromosome, position);

                // --------------------- If genotypes is empty, then this position is missing within the index... Skip ahead.
                if ! fst_reader.has_genotypes() {
                    debug!("Missing coordinate in fst index: [{chromosome: <2} {position: >9}]");
                    continue 'coordinate
                }

                // --------------------- Skip line if allele frequency is < maf and keep in memory for filtration..
                let pop_af = fst_reader.get_pop_allele_frequency(&self.pedigree_pop)?;
                if pop_af < maf || pop_af > (1.0-maf) { 
                    debug!("skip allele at [{chromosome: <2} {position: >9}]: pop_af: {pop_af:<8.5} --maf: {maf}");
                    positions_to_delete.entry(key.to_owned()).or_insert_with(|| Vec::new()).push(pairwise_diff.coordinate);
                    continue 'coordinate
                }
                
                let pileup_error_probs = pairwise_diff.error_probs();
                // --------------------- Parse genotype fields and start updating dynamic simulations.
                self.update_pedigrees(&fst_reader, chromosome, position, &key, &pileup_error_probs)?;
                //true
            }
        }

        // --------------------- Filter out unwanted alleles. 
        Self::filter(&mut positions_to_delete, comparisons);
        Ok(())
    }

    fn filter(positions_to_delete: &mut HashMap<String, Vec<Coordinate>>, comparisons: &mut Comparisons) {
        info!("Filtering out unwanted alleles from comparisons.");
        for comparison in comparisons.iter_mut() {
            let key = comparison.get_pair();
            let pre_filtered_n = comparison.positions.len();
            comparison.positions.retain(|pwd| ! positions_to_delete.entry(key.to_owned()).or_default().contains(&&pwd.coordinate));
            info!("- {key: <20} filtered-out SNPs: {}", pre_filtered_n - comparison.positions.len());
        }
    }

    pub fn pedigree_simulations_vcf(&mut self, comparisons: &mut Comparisons, input_vcf_path: &Path, maf: f64, threads: usize) -> Result<(), Box<dyn Error>> {
        // --------------------- Read VCF File line by line.
        let mut i = 0;
        let mut vcf_reader = VCFReader::new(input_vcf_path, threads)?;

        // Keep a record of positions that should get filtered out after maf < treshold.
        let mut positions_to_delete = HashMap::new();

        'line: while vcf_reader.has_data_left()? {

            // Get current chromosome and position.
            let (chromosome, position) = vcf_reader.parse_coordinate()?; // 1
            let coordinate = Coordinate::new(chromosome, position);


            // --------------------- Print progress in increments of 50 000 lines.
            if i % 50_000 == 0 { info!(" {i: >9}: [{chromosome: <2} {position: >9}]");}
            i+=1;

            // --------------------- Loop across comparisons. 
            //let mut genotypes_filled: bool = false;
            'comparison: for comparison in comparisons.iter(){
                let key = comparison.get_pair();
                // --------------------- Skip if the current position is not a valid candidate.
                let relevant_position = comparison.positions.get(&coordinate);
                match relevant_position {
                    None                      => { continue 'comparison },
                    Some(pairwise_diff) => {
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
                                    trace!("Skip allele at [{chromosome: <2} {position: >9}]: pop_af: {pop_af:<8.5} --maf: {maf}");
                                    positions_to_delete.entry(key.to_owned()).or_insert_with(|| Vec::new()).push(pairwise_diff.coordinate);
                                    vcf_reader.next_line()?;
                                    continue 'line
                            }

                            // Parse genotypes. 
                            vcf_reader.fill_genotypes()?;
                        }
                        // Parse genotype fields and start updating dynamic simulations.
                        let pileup_error_probs = pairwise_diff.error_probs();
                        self.update_pedigrees(&vcf_reader, chromosome, position, &comparison.get_pair().to_owned(), &pileup_error_probs)?;          
                    }
                }
            }
            // Reset line if we never parsed genotypes.
            vcf_reader.next_line()?;
        }

        // --------------------- Filter out unwanted alleles. 
        Self::filter(&mut positions_to_delete, comparisons);
        Ok(())        
    }

    pub fn write_simulations(&self, comparisons: &Comparisons, output_files: &HashMap<String, String>) -> Result<(), Box<dyn Error>> {
        // --------------------- Print pedigree simulation results.
        for comparison in comparisons.iter() {
            let comparison_label = comparison.get_pair();
            let pedigree_vec = self.pedigrees.get(&comparison_label).unwrap();

            let mut writer = pwd_from_stdin::io::Writer::new(Some(output_files[&comparison_label].clone()))?;
            writer.write_iter(vec![&pedigree_vec])?;

            info!("\n--------------------- {comparison_label}\n{}", pedigree_vec);
        }
        Ok(())
    }

    pub fn compute_results(&self, comparisons: &Comparisons, output_file: &String) -> Result<(), Box<dyn Error>> {
        let mut simulations_results = Vec::with_capacity(comparisons.len());
        let mut writer = pwd_from_stdin::io::Writer::new(Some(output_file.clone()))?;

        // Print header and write to result file
        let simulation_header = format!("{: <20} - {: <20} - {: <10} - {: <10} - {: <10}",
            "Pair_name", "Most_Likely_rel.", "Corr.avg-PWD", "Sim.avg-PWD", "Min_Z-Score"
        );
        println!("{simulation_header}");
        simulations_results.push(simulation_header);


        for comparison in comparisons.iter() {
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

            let simulation_result = format!("{pair_name: <20} - {most_likely_rel: <20} - {observed_avg_pwd: >12.6} - {most_likely_avg_pwd: >12.6} - {min_z_score: <12.6}");
            println!("{simulation_result}");
            simulations_results.push(simulation_result);
        }


        writer.write_iter(simulations_results)?;

        Ok(())
    }
}