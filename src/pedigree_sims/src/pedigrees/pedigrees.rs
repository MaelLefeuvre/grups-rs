use std::{
    collections::HashMap,
    path::{PathBuf, Path},
    error::Error
};

use super::{
    PedigreeReps,
    pedigree::pedigree_params::ParamRateGenerator,
};

use crate::io::{
    self,
    genotype_reader::GenotypeReader,
    vcf::{
        SampleTag,
        reader::{VCFReader, VCFPanelReader},
    }, 
};

use genome::GeneticMap;

use pwd_from_stdin::{
    self,
    comparisons::{
        Comparisons,
        pwd::Coordinate,
    },
};

use rand::{Rng, prelude::ThreadRng};
use log::{info, trace, debug};


/// A Pedigree simulator, aggregating all pileup-comparison pedigree simulation replicates.
/// # Fields
/// - `pedigrees`         : HashMap of pedigree simulation replicates.
///                         Key = pileup comparison label | value = pedigree simulation replicates
/// - `pedigree_pop`      : (super)-population id used for the pedigree simulations.
/// - `genetic_map`       : Genetic recombination map interval tree.
/// - `previous_positions`: Hashmap, tracking the coordinate of each previously typed SNP coordinate, for a given pileup comparison
///                         Key = pileup comparison label | value = position of the previous SNP coordinate
/// - `rng`               : random number generator.
/// 
/// # TODO:
///  - `previous_position` is not fit to handle multiple chromosomes at this state... This could cause bugs if one of the input files
///    contains multiple chromosomes... 
pub struct Pedigrees {
    pedigrees: HashMap<String, PedigreeReps>,
    pedigree_pop: String,
    genetic_map: GeneticMap,
    previous_positions: HashMap<String, u32>,
    rng: ThreadRng,
}

/// @TODO! : Right now we're reading the pedigree definition file for each replicate. Which isn't very efficient.
///          This is because simply using clone() on a pedigree_template will merely clone() the Rc<RefCell>.
///          ==> This is NOT what we want. What we want is to create new instances of Rc<RefCell> from these values.
impl Pedigrees {
    /// Initialize a pedigree simulator.
    /// # Arguments:
    /// - `pedigree_pop`: (super-)population id used for the pedigree simulation replicates.
    /// - `comparisons` : pileup Comparisons of our real samples.
    /// - `recomb_dir`  : path leading to the directory containing recombination maps.
    /// 
    /// # Errors:
    /// - returns an error upon failing to parse `self.genetic_map`
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

    /// Iterate upon the pileup Comparisons of our real samples, and initialize `n` pedigree simulation replicates for them.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `panel`         : input samples panel definition file.
    /// - `reps`          : user-defined number of pedigree simulation replicates.
    /// - `pedigree_path` : path leading to the pedigree definition file
    /// - contam_pop`     : user-defined name of the (super-)population-id requested to simulate modern human contamination
    /// - `contam_num_ind`: vector of user-requested contaminating individuals `contam_num_ind[i] is tied to pileup `sample[i]`
    pub fn populate(&mut self, comparisons: &Comparisons,  panel: &VCFPanelReader, reps: u32, pedigree_path: &Path, contam_pop: &Vec<String>, contam_num_ind: &Vec<usize>) -> Result<(), Box<dyn Error>> {

        // ---- Randomly sample contaminant SampleTags, according to the requested population and number of individuals.
        let samples_contam_tags: Vec<Vec<SampleTag>> = panel.fetch_contaminants(contam_pop, contam_num_ind)?;

        // ---- Iterate on each pileup comparison, and initialize pedigree simulation replicates for them.
        for comparison in comparisons.iter() {
            let comparison_label = comparison.get_pair();
            let pair_indices = comparison.get_pair_indices();

            let mut pedigree_reps = PedigreeReps::with_capacity(reps as usize);

            // ---- Assign contaminating individuals to the pedigree vector.
            pedigree_reps.set_contaminants(&samples_contam_tags, pair_indices);
            debug!("Contaminant set for {comparison_label}:\n{}", pedigree_reps.contaminants.as_ref().unwrap());

            // ---- Initialize pedigree replicates.
            pedigree_reps.populate(pedigree_path, &self.pedigree_pop, panel)?;
            self.pedigrees.insert(comparison_label.to_owned(), pedigree_reps);
        }
        Ok(())
    }

    /// Iterate upon the pileup comparisons of our real samples and define simulation parameters for each pedigree 
    /// simulation replicate, using the user input.
    /// Arguments:
    /// - `comparisons`           : pileup Comparisons of our real samples.
    /// - `snp_downsampling_rate` : user-defined probability of ignoring an snp coordinate during simulations
    /// - `af_downsampling_rate`  : user-defined probability of performing allele fixation during simulations
    /// - `seq_error_rates`       : user-defined probabilities of simulating a sequencing error (constant values and/or ranges)
    ///                             `seq_error_rate[i]` is tied to pileup `sample[i]`
    /// - `seq_error_rates`       : user-defined probabilities of simulating a human contamination (constant values and/or ranges)
    ///                             `contam_rate[i]` is tied to pileup `sample[i]`
    /// 
    /// # Errors
    /// - if `self.pedigrees` does contain a given pileup comparison label.
    pub fn set_params(&mut self, comparisons: &Comparisons, snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: &Option<Vec<Vec<f64>>>, contam_rate: &Vec<Vec<f64>>) -> Result<(), String>{

        // ---- Iterate upon pileup comparisons and assign parameters for each pedigree simulation replicate.
        for comparison in comparisons.iter() {
            let pair_indices = comparison.get_pair_indices();
            let pair_label = comparison.get_pair();
            let pedigree_reps = self.pedigrees.get_mut(&pair_label).ok_or(format!("Cannot access pair {pair_label}"))?;

            // ---- Instantiate a sequencing error `ParamRateGenerator` if the user specified sequencing error rates.
            //      If the user did not provide any, assign `None` -> the phred-scores of the pileup will then be used to compute the seq-error probability
            let mut seq_error_rate_gen = match seq_error_rate {
                Some(seq_errors_vec) => Some(ParamRateGenerator::from_user_input(seq_errors_vec, pair_indices)),
                None => None,
            };
            // ---- Instantiate a contamination `ParamRateGenerator`
            let mut contam_rate_gen = ParamRateGenerator::from_user_input(contam_rate, pair_indices);

            // ---- Iterate on each pedigree replicate, and use our two `ParamRateGenerator` to assign random and/or
            //      constant values. (depending on the user-input)
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

    /// Iterate upon all the pedigree replicates of a given pileup comparison, and update our simulated avg_pwd.
    /// Arguments:
    /// - `reader`            : a `GenotypeReader` trait object. Either `VCFReader` or `FSTReader`
    /// - `chromosome`        : chromosome name
    /// - `position`          : SNP position within the chromosome.
    /// - `comparison_label`  : pileup comparison label (e.g. 'Ind1-Ind2')
    /// - `pileup_error_probs`: local pileup sequencing error probabilities. Extracted from phred-scores.
    /// 
    /// # Errors
    /// - if any pedigree's `self.params` field is set to `None`
    /// - if any pedigree vector's contaminant is set to `None`
    /// # Panics:
    ///  - if `comparison_label` does not match any key within `self.pedigrees`
    ///  - whenever a contaminating individual carries multi-allelic alternate allele (i.e. alt allele is > 1)
    fn update_pedigrees(&mut self, reader: &dyn GenotypeReader, chromosome:u8, position: u32, comparison_label: &String, pileup_error_probs: &[f64;2]) -> Result<(), Box<dyn Error>>{

        // ---- Compute the interval between current and previous position, search trough the genetic map interval tree,
        //      and compute the probability of recombination.
        let previous_position    : u32 = self.previous_positions[comparison_label];
        let interval_prob_recomb : f64 = self.genetic_map.compute_recombination_prob(chromosome, previous_position, position);
        //trace!("SNP candidate for {comparison_label} - [{chromosome:<2} {position:>9}] - pop_af: {pop_af:?} - cont_af: {cont_af:?} - recomb_prop: {interval_prob_recomb:<8.6}");


        let pedigree_vec = self.pedigrees.get_mut(comparison_label).unwrap();

        // -------------------- Get the contaminating population allele frequency
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

    /// Perform pedigree simulations using an `FSTReader`.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `input_fst_path`: path leading to the target FST-index file, containing genotype information for founder individuals. (`.fst`).
    /// - `maf`           : user-defined minor-allele-frequency treshold.
    /// 
    /// # Panics:
    /// - when failing to convert `input_fst_path` to a string slice.
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
            }
        }

        // --------------------- Filter out unwanted alleles. 
        Self::filter_pileup_positions(&mut positions_to_delete, comparisons);
        Ok(())
    }

    /// Filter out a subset of target coordinates from our pileup comparisons.
    /// This is mainly used to remove positions that are under the user-defined `--maf` treshold, and thus recompute a
    /// 'corrected' observed avg. pairwise difference.
    /// # Arguments:
    /// - `positions_to_delete`: HashMap of coordinates. Key = comparison label | Value = Vec<Coordinate> to delete.
    /// - `comparisons`        : pileup Comparisons of our real samples.
    fn filter_pileup_positions(positions_to_delete: &mut HashMap<String, Vec<Coordinate>>, comparisons: &mut Comparisons) {
        info!("Filtering out unwanted alleles from comparisons.");
        for comparison in comparisons.iter_mut() {
            let key = comparison.get_pair();
            let pre_filtered_n = comparison.positions.len();
            comparison.positions.retain(|pwd| ! positions_to_delete.entry(key.to_owned()).or_default().contains(&&pwd.coordinate));
            info!("- {key: <20} filtered-out SNPs: {}", pre_filtered_n - comparison.positions.len());
        }
    }

    /// Perform pedigree simulations using a `VCFReader`.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `input_vcf_path`: path leading to the target VCF file, containing genotype information for founder individuals. (`.vcf` or `.vcf.gz`).
    /// - `maf`           : user-defined minor-allele-frequency treshold.
    /// - `threads`       : user-requested number of decompression threads for our VCFReader
    ///                    (this is only relevant when reading BGZF-compressed `.vcf.gz` files)
    /// 
    /// # Panics:
    /// - when failing to convert `input_vcf_path` to a string slice.
    pub fn pedigree_simulations_vcf(&mut self, comparisons: &mut Comparisons, input_vcf_path: &Path, maf: f64, threads: usize) -> Result<(), Box<dyn Error>> {
        // --------------------- Read VCF File line by line.
        let mut i = 0;
        let mut vcf_reader = VCFReader::new(input_vcf_path, threads)?;

        // --------------------- Keep a record of positions that should get filtered out
        //                       after maf < treshold.
        let mut positions_to_delete = HashMap::new();

        'line: while vcf_reader.has_data_left()? {

            // --------------------- Get current chromosome and position.
            let (chromosome, position) = vcf_reader.parse_coordinate()?; // 1
            let coordinate = Coordinate::new(chromosome, position);


            // --------------------- Print progress in increments of 50 000 lines.
            if i % 50_000 == 0 { info!(" {i: >9}: [{chromosome: <2} {position: >9}]");}
            i+=1;

            // --------------------- Loop across comparisons. 
            'comparison: for comparison in comparisons.iter(){
                let key = comparison.get_pair();
                // --------------------- Skip if the current position is not a valid candidate.
                let relevant_position = comparison.positions.get(&coordinate);
                match relevant_position {
                    None                      => { continue 'comparison },
                    Some(pairwise_diff) => {
                        if ! vcf_reader.genotypes_filled {
                            // ---- Go to INFO field.
                            vcf_reader.parse_info_field()?;

                            // ---- Check if this coordinate is a biallelic SNP and skip line if not.
                            if vcf_reader.is_multiallelic() || ! vcf_reader.is_snp() {
                                vcf_reader.next_line()?;
                                continue 'line
                            }

                            // ---- Extract population allele frequency.
                            let pop_af = vcf_reader.get_pop_allele_frequency(&self.pedigree_pop)?;

                            // ---- Skip line if allele frequency is < maf
                            //      [WARN]: Note that {POP}_AF entries in the INFO field relate to the REF allele frequency, not MAF.
                            //              ==> for maf = 0.05, we must also filter out alleles with AF > 0.95.
                            if pop_af < maf || pop_af > (1.0-maf) {
                                    trace!("Skip allele at [{chromosome: <2} {position: >9}]: pop_af: {pop_af:<8.5} --maf: {maf}");
                                    positions_to_delete.entry(key.to_owned()).or_insert_with(|| Vec::new()).push(pairwise_diff.coordinate);
                                    vcf_reader.next_line()?;
                                    continue 'line
                            }

                            // ---- Parse genotypes. 
                            vcf_reader.fill_genotypes()?;
                        }
                        // ---- Parse genotype fields and start updating dynamic simulations.
                        let pileup_error_probs = pairwise_diff.error_probs();
                        self.update_pedigrees(&vcf_reader, chromosome, position, &comparison.get_pair().to_owned(), &pileup_error_probs)?;          
                    }
                }
            }
            // ---- Reset line if we never parsed genotypes.
            vcf_reader.next_line()?;
        }

        // --------------------- Filter out unwanted alleles. 
        Self::filter_pileup_positions(&mut positions_to_delete, comparisons);
        Ok(())        
    }

    /// Log pedigree simulation results and write them into a set of output files.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `output_files`  : HashMap of target output files.
    ///                     Key = comparison label | value = target output path.
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

    /// Recompute a corrected average-PWD of our real samples, and estimate the most likely relationship
    /// using our simulation results
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `output_files`  : target output file where results are written.
    pub fn compute_results(&self, comparisons: &Comparisons, output_file: &String) -> Result<(), Box<dyn Error>> {
        // ---- Resource acquisition
        let mut simulations_results = Vec::with_capacity(comparisons.len());
        let mut writer = pwd_from_stdin::io::Writer::new(Some(output_file.clone()))?;

        // ---- Print header and write to output_file
        let simulation_header = format!("{: <20} - {: <20} - {: <10} - {: <10} - {: <10}",
            "Pair_name", "Most_Likely_rel.", "Corr.avg-PWD", "Sim.avg-PWD", "Min_Z-Score"
        );
        println!("{simulation_header}");
        simulations_results.push(simulation_header);

        // ---- loop across our pileup comparisons and assign a most-likely relationship, using our simulations.
        for comparison in comparisons.iter() {
            let comparison_label    = comparison.get_pair();
            let pedigree_vec = self.pedigrees.get(&comparison_label).unwrap();
    
            // ---- Aggregate the sum of avg. PWD for each relatedness scenario.
            let sum_simulated_pwds = pedigree_vec.compute_sum_simulated_pwds();
    
            // ---- Select the scenario having the least amount of Z-score with our observed avg.PWD
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
    
            // ---- Correct with f64::NaN in case there were not enough SNPs to compute a simulated avg.pwd.
            let min_z_score = if min_z_score == f64::MAX {f64::NAN} else {min_z_score};
            let pair_name = comparison.get_pair();

            // ---- Preformat and log result to console.
            let simulation_result = format!("{pair_name: <20} - {most_likely_rel: <20} - {observed_avg_pwd: >12.6} - {most_likely_avg_pwd: >12.6} - {min_z_score: <12.6}");
            println!("{simulation_result}");
            simulations_results.push(simulation_result);
        }

        // ---- Write simulation results to file.
        writer.write_iter(simulations_results)?;

        Ok(())
    }
}