use std::{collections::HashMap, path::Path, mem::ManuallyDrop};


use grups_io::{
    write::GenericWriter,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader, VCFReader, FSTReader, fst::SetRead},
    read::PanelReader,
};

use located_error::prelude::*;
use genome::{coordinate::{Coordinate, Position}, GeneticMap};
use parser::RelAssignMethod;
use pwd_from_stdin::comparisons::{Comparisons as PileupComparisons, Comparison};

use fastrand;
use ahash::AHashMap;
use indexmap::IndexMap;
use log::{info, trace, debug, warn};

mod pedigree_reps;
pub use pedigree_reps::PedigreeReps;
mod pedigree;
use pedigree::{Pedigree, pedparam::ParamRateGenerator};
use pedigree::Contaminant;

mod error;
use error::PedigreeError;

mod pedigree_sim_builder;

pub mod constants;
use constants::{
    REPLICATE_ID_FORMAT_LEN,
    COMPARISON_LABEL_FORMAT_LEN,
    IND_LABEL_FORMAT_LEN,
    IND_TAG_FORMAT_LEN,
    PWD_FORMAT_LEN,
    OVERLAP_FORMAT_LEN,
    AVG_PWD_FORMAT_LEN,
    SEX_FORMAT_LEN,
    FLOAT_FORMAT_PRECISION
};

use crate::svm::LibSvmBuilder;


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
    pedigrees         : HashMap<String, PedigreeReps>,
    pedigree_pop      : String,
    genetic_map       : GeneticMap,
    previous_positions: HashMap<String, Position>,
    rng               : fastrand::Rng,
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
    pub fn initialize(pedigree_pop: &str, comparisons: &PileupComparisons, recomb_dir: impl AsRef<Path>) -> Result<Self> {
        // Generate pedigree replicates for each pwd_from_stdin::Comparison.
        let pedigrees = HashMap::new();


        // --------------------- Parse input recombination maps.
        //info!("Parsing genetic maps in {}", recomb_dir);
        let genetic_map = GeneticMap::from_dir(recomb_dir).loc("While attempting to initialize Pedigrees")?;

        // --------------------- For each comparison, keep a record of the previously typed SNP's position.
        let mut previous_positions = HashMap::new();
        for comparison in comparisons.iter() {
            let key = comparison.get_pair();
            previous_positions.insert(key.to_owned(), Position(0));
        }
        // --------------------- Initialize RNG
        let rng = fastrand::Rng::new();
        let pedigree_pop = pedigree_pop.to_string();
        Ok(Pedigrees{pedigrees, pedigree_pop, previous_positions, genetic_map, rng})
    }

    pub fn set_founder_tags(&mut self, panel: &PanelReader) -> Result<()> {
        self.pedigrees.iter_mut().try_for_each(|(label, ped_rep)|{
            ped_rep.set_founder_tags(panel, &self.pedigree_pop)
                .with_loc(||format!("While attempting to set founder tags in pedigree vector of comparison '{label}'"))
        })
    }

    pub fn assign_offspring_strands(&mut self) -> Result<()> {
        self.pedigrees.iter_mut().try_for_each(|(label, ped_rep)|{
            ped_rep.assign_offspring_strands()
                .with_loc(||format!("While attempting to randomly assign offspring strands in pedigree vector of comparison '{label}'"))

        })
    }

    /// Iterate upon the pileup Comparisons of our real samples, and initialize `n` pedigree simulation replicates for them.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `panel`         : input samples panel definition file.
    /// - `reps`          : user-defined number of pedigree simulation replicates.
    /// - `pedigree_path` : path leading to the pedigree definition file
    /// - contam_pop`     : user-defined name of the (super-)population-id requested to simulate modern human contamination
    /// - `contam_num_ind`: vector of user-requested contaminating individuals `contam_num_ind[i] is tied to pileup `sample[i]`
    pub fn populate(&mut self, comparisons: &PileupComparisons,  panel: &PanelReader, reps: u32, pedigree_path: &Path, contam_pop: &[String], contam_num_ind: &Vec<usize>) -> Result<()> {

        // ---- Randomly sample contaminant SampleTags, according to the requested population and number of individuals.
        let samples_contam_tags: Vec<Vec<SampleTag>> = panel.fetch_contaminants(contam_pop, contam_num_ind)
            .with_loc(|| format!("Failed to fetch contaminating individuals from the {contam_pop:?} population"))?;

        // ---- Iterate on each pileup comparison, and initialize pedigree simulation replicates for them.
        for comparison in comparisons.iter() {
            let comparison_label = comparison.get_pair();
            let pair_indices = comparison.get_pair_indices();

            let mut pedigree_reps = PedigreeReps::with_capacity(reps as usize);

            // ---- Assign contaminating individuals to the pedigree vector.
            pedigree_reps.set_contaminants(&samples_contam_tags, pair_indices);

            // ---- Debug print.
            debug!("Contaminant set for {comparison_label}:\n{}", 
                match pedigree_reps.contaminants.as_ref() {
                    Some(cont) => cont.to_string(),
                    None => String::from("None")
                }
            );

            // ---- Initialize pedigree replicates.
            pedigree_reps.populate(pedigree_path)
                .map_err(PedigreeError::PopulateError)
                .with_loc(|| format!("While attempting to populate a vector of pedigrees for {}", comparison.get_pair()))?;
            self.pedigrees.insert(comparison_label.to_owned(), pedigree_reps);
        }
        Ok(())
    }

    pub fn all_sex_assigned(&self) -> bool {
        self.pedigrees.values().all(|ped_rep| ped_rep.all_sex_assigned())
    }

    pub fn assign_random_sex(&mut self) -> Result<()> {
        self.pedigrees.iter_mut().try_for_each(|(label, ped_rep)| {
            ped_rep.assign_random_sex().with_loc(|| format!("While assigning sexes in pedigree replicate vector of comparison {label}"))
        })
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
    pub fn set_params(&mut self, 
        comparisons          : &PileupComparisons, 
        snp_downsampling_rate: f64,
        af_downsampling_rate : f64,
        seq_error_rate       : &Option<Vec<Vec<f64>>>,
        contam_rate          : &[Vec<f64>]
    ) -> Result<()>{
        let loc_msg = "While attempting to set user input parameters";
        // ---- Iterate upon pileup comparisons and assign parameters for each pedigree simulation replicate.
        for comparison in comparisons.iter() {
            use PedigreeError::MissingPedVec;
            let pair_indices  = comparison.get_pair_indices();
            let pair_label    = comparison.get_pair();
            let pedigree_reps = self.pedigrees.get_mut(&pair_label).ok_or_else(|| MissingPedVec(pair_label.to_string())).loc(loc_msg)?;

            // ---- Instantiate a sequencing error `ParamRateGenerator` if the user specified sequencing error rates.
            //      If the user did not provide any, assign `None` -> the phred-scores of the pileup will then be used to compute the seq-error probability
            let mut seq_error_rate_gen = seq_error_rate.as_ref()
                .map(|seq_errors_vec| ParamRateGenerator::from_user_input(seq_errors_vec, pair_indices));

            // ---- Instantiate a contamination `ParamRateGenerator`
            let mut contam_rate_gen = ParamRateGenerator::from_user_input(contam_rate, pair_indices);

            match seq_error_rate_gen {
                Some(ref seq_err) => debug!("  - seq-error-rate for pair {}: {seq_err}", comparison.get_pair()),
                None              => debug!("  - seq-error-rate for pair {}: pileup", comparison.get_pair())
            };
            debug!("  - contam-rate    for pair {}: {contam_rate_gen}", comparison.get_pair());

            // ---- Iterate on each pedigree replicate, and use our two `ParamRateGenerator` to assign random and/or
            //      constant values. (depending on the user-input)
            pedigree_reps.iter_mut()
                .for_each(|ped| {
                    let seq_error_rate = seq_error_rate_gen.as_mut()
                    .map(|generator| generator.gen_random_values());
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
    // #[inline(always)]
    fn update_pedigrees(&mut self, reader: &dyn GenotypeReader, coordinate: &Coordinate, comparison_label: &String, pileup_error_probs: &[f64;2]) -> Result<()>{
        use PedigreeError::{InvalidCoordinate, MissingContaminant, MissingPedVec};
        // ---- Compute the interval between current and previous position, search trough the genetic map interval tree,
        //      and compute the probability of recombination.
        let previous_position    = self.previous_positions.get_mut(comparison_label).with_loc(||InvalidCoordinate)?;
        let interval_prob_recomb = self.genetic_map.compute_recombination_prob(coordinate, *previous_position);
        trace!("SNP candidate for {comparison_label} - {coordinate} - recomb_prob: {interval_prob_recomb:<8.6}");

        // ---- Extract the vector of pedigrees that'll get updated.
        //      @TODO: This error handling is performed multiple times. Stay DRY & wrap this in a public method.
        let pedigree_vec = self.pedigrees.get_mut(comparison_label)
            .ok_or_else(|| MissingPedVec(comparison_label.to_string()))
            .loc("While attempting to update pedigrees")?;
        // -------------------- Get the contaminating population allele frequency
        let cont_af = pedigree_vec.contaminants.as_ref().with_loc(||MissingContaminant)?.compute_local_cont_af(reader)?;

        let xchr_mode = coordinate.chromosome.0 == b'X';
        'pedigree: for (i, pedigree) in pedigree_vec.iter_mut().enumerate() {
            // --------------------- Perform SNP downsampling if necessary
            if self.rng.f64() < pedigree.get_params()?.snp_downsampling_rate {continue 'pedigree}

            // --------------------- Update founder alleles. Perform allele frequency downsampling if necessary.
            pedigree.update_founder_alleles(reader, &mut self.rng)?;

            // --------------------- Compute offspring genomes
            pedigree.compute_offspring_alleles(interval_prob_recomb, i, &mut self.rng, xchr_mode)?;

            // --------------------- Compare genomes.
            pedigree.compare_alleles(cont_af, pileup_error_probs, &mut self.rng)?;
            // --------------------- Clear genotypes before the next line!
            pedigree.clear_alleles();
        }
        // Keep track of the previous position:
        *previous_position = coordinate.position;
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
    pub fn pedigree_simulations_fst<T: SetRead<T> + AsRef<[u8]>>(&mut self, comparisons: &mut PileupComparisons, input_fst_path: &Path, maf: f32) -> Result<()> {
        let loc_msg = "While performing pedigree simulations";
        // ----  Initialize a new FSTReader.
        let mut fst_reader: FSTReader::<T> = FSTReader::new(&input_fst_path.to_string_lossy()).loc(loc_msg)?;

        // ---- Get a list of which chromosomes are found within the FST index. 
        let contained_chromosomes = fst_reader.find_chromosomes().loc(loc_msg)?;
        debug!("Contained chromosomes: {contained_chromosomes:?}");
        
        // ---- Extract the minimum and maximum values.
        let max = *contained_chromosomes.last().context("Failed to retrieve chromosome range").loc(loc_msg)?;
        let min = contained_chromosomes[0];

        // ---- Convert these values as a range of valid coordinates.
        let start = Coordinate::new(min, std::u32::MIN);
        let stop  = Coordinate::new(max, std::u32::MAX);

        // Keep a record of positions that should get filtered out after maf < treshold.
        let mut positions_to_delete: AHashMap<String, Vec<Coordinate>> = AHashMap::new();
        'comparison: for comparison in comparisons.iter() {
            let mut missing_snps: u32 = 0; //MODIFIED
            info!("Performing simulations for : {}", comparison.get_pair());

            // ---- Extract positions matching the FST index chromosome(s)
            // ---- NOTE: Here we're using range(). But what we'd ideally like is drain_filter() . [Nightly only]
            let relevant_positions = comparison.positions.range(start..stop);
            //
            // Get the number of typed positions for these chromosomes (this .clone() is quite expensive for just a fancy progress estimation...)
            let n = relevant_positions.clone().count();

            // Loop along positions.
            let key = comparison.get_pair();
            'coordinate: for (i, pairwise_diff) in relevant_positions.enumerate() {                
                let coordinate = pairwise_diff.coordinate;
                // --------------------- Print progress in increments of 10%
                if (i % ((n/10)+1)) == 0 {
                    let percent = (i as f32 / (n as f32).floor()) * 100.0 ;
                    info!("{percent: >5.1}% : {coordinate}");
                }

                // --------------------- Search through the FST index for the genotypes and pop frequencies at this coordinate.
                fst_reader.clear_buffers();
                fst_reader.search_coordinate_genotypes(&coordinate);
                fst_reader.search_coordinate_frequencies(&coordinate);

                // --------------------- If genotypes is empty, then this position is missing within the index... Skip ahead.
                if ! fst_reader.has_genotypes() {
                    trace!("Missing coordinate in fst index: {coordinate}");
                    missing_snps += 1;                    // MODIFIED
                    continue 'coordinate
                }

                // --------------------- Skip line if allele frequency is < maf and keep in memory for filtration..
                let pop_af = fst_reader.get_pop_allele_frequency(&self.pedigree_pop)?;
                if pop_af < maf || pop_af > (1.0-maf) { 
                    trace!("skip allele at {coordinate}: pop_af: {pop_af:<8.5} --maf: {maf}");
                    positions_to_delete.entry(key.to_owned()).or_default().push(pairwise_diff.coordinate);
                    continue 'coordinate
                }
                
                let pileup_error_probs = pairwise_diff.error_probs();
                // --------------------- Parse genotype fields and start updating dynamic simulations.
                self.update_pedigrees(&fst_reader, &coordinate,  &key, &pileup_error_probs)?;
            }

            // ---- Count missing SNPs as non-informative positions. i.e. an overlap.
            if missing_snps > 0 {
                info!("{missing_snps} missing SNPs within the simulation dataset. Counting those as non-informative overlap.");
                self.pedigrees.get_mut(&comparison.get_pair()).unwrap().add_non_informative_snps(missing_snps);
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
    fn filter_pileup_positions(positions_to_delete: &mut AHashMap<String, Vec<Coordinate>>, comparisons: &mut PileupComparisons) {
        info!("Filtering out unwanted alleles from comparisons.");
        for comparison in comparisons.iter_mut() {
            let key = comparison.get_pair();
            let pre_filtered_n = comparison.positions.len();
            comparison.positions.retain(|pwd| ! positions_to_delete.entry(key.to_owned()).or_default().contains(&pwd.coordinate));
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
    pub fn pedigree_simulations_vcf(&mut self, comparisons: &mut PileupComparisons, input_vcf_path: &Path, maf: f32, threads: usize) -> Result<()> {
        let loc_file  = || format!("While attempting to perform pedigree simulations on {}", input_vcf_path.display());
        let loc_coord = {|c: &Coordinate| format!("While parsing coordinate: {c}")};
        // --------------------- Read VCF File line by line.
        let mut i = 0;
        let mut vcf_reader = VCFReader::new(input_vcf_path, threads)
            .with_loc(loc_file)?;


        // --------------------- Keep a record of positions that should get filtered out
        //                       after maf < treshold.
        let mut positions_to_delete: AHashMap<String, Vec<Coordinate>> = AHashMap::new();
        'line: while vcf_reader.has_data_left().with_loc(loc_file)? {

            // --------------------- Get current chromosome and position.
            let coordinate = vcf_reader.parse_coordinate().with_loc(loc_file)?;

            // --------------------- Print progress in increments of 50 000 lines.
            if i % 50_000 == 0 { info!(" {i: >9}: {coordinate}");}
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
                            vcf_reader.parse_info_field().with_loc(|| loc_coord(&coordinate))?;

                            // ---- Check if this coordinate is a biallelic SNP and skip line if not.
                            let not_an_snp = !vcf_reader.is_snp().with_loc(|| loc_coord(&coordinate))?;
                            if vcf_reader.is_multiallelic() || not_an_snp {
                                vcf_reader.next_line().with_loc(|| loc_coord(&coordinate))?;
                                continue 'line
                            }

                            // ---- Extract population allele frequency.
                            let pop_af = vcf_reader.get_pop_allele_frequency(&self.pedigree_pop)
                                .with_loc(|| loc_coord(&coordinate))?;

                            // ---- Skip line if allele frequency is < maf
                            //      [WARN]: Note that {POP}_AF entries in the INFO field relate to the REF allele frequency, not MAF.
                            //              ==> for maf = 0.05, we must also filter out alleles with AF > 0.95.
                            if pop_af < maf || pop_af > (1.0-maf) {
                                    trace!("Skip allele at {coordinate}: pop_af: {pop_af:<8.5} --maf: {maf}");
                                    positions_to_delete.entry(key.to_owned()).or_default().push(pairwise_diff.coordinate);
                                    vcf_reader.next_line().with_loc(|| loc_coord(&coordinate))?;
                                    continue 'line
                            }

                            // ---- Parse genotypes. 
                            vcf_reader.fill_genotypes().with_loc(|| loc_coord(&coordinate))?;
                        }
                        // ---- Parse genotype fields and start updating dynamic simulations.
                        let pileup_error_probs = pairwise_diff.error_probs();
                        self.update_pedigrees(&vcf_reader, &coordinate, &comparison.get_pair().to_owned(), &pileup_error_probs)
                            .with_loc(|| loc_coord(&coordinate))?;          
                    }
                }
            }
            // ---- Reset line if we never parsed genotypes.
            vcf_reader.next_line().with_loc(|| loc_coord(&coordinate))?;
        }

        // --------------------- Filter out unwanted alleles. 
        Self::filter_pileup_positions(&mut positions_to_delete, comparisons);
        Ok(())        
    }

    fn get_pedigree_vec(&self, comparison_label: &str) -> Result<&PedigreeReps, PedigreeError> {
        use PedigreeError::MissingPedVec;
        self.pedigrees.get(comparison_label).ok_or_else(|| MissingPedVec(comparison_label.to_string()))
    }

    /// Log pedigree simulation results and write them into a set of output files.
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `output_files`  : HashMap of target output files.
    ///                     Key = comparison label | value = target output path.
    pub fn write_simulations(&self, comparisons: &PileupComparisons, output_files: &HashMap<String, String>) -> Result<()> {
        let loc_msg = "While attempting to write simulation results";

        // --------------------- Print pedigree simulation results.
        trace!("Printing pairwise simulation results:");
        let simulations_header = format!(
            "{: <REPLICATE_ID_FORMAT_LEN$} - \
            {: <COMPARISON_LABEL_FORMAT_LEN$} - \
            {: <IND_LABEL_FORMAT_LEN$} - \
            {: <IND_LABEL_FORMAT_LEN$} - \
            {: <IND_TAG_FORMAT_LEN$} - \
            {: <IND_TAG_FORMAT_LEN$} - \
            {: >PWD_FORMAT_LEN$} - \
            {: >OVERLAP_FORMAT_LEN$} - \
            {: <AVG_PWD_FORMAT_LEN$.FLOAT_FORMAT_PRECISION$} - \
            {: <SEX_FORMAT_LEN$} - \
            {: <SEX_FORMAT_LEN$}",
            "replicate", "label", "parent0", "parent1", "parent0.id", "parent1.id", "pwd", "overlap", "avg", "parent0.sex", "parent1.sex"
        );
        for comparison in comparisons.iter() {
            let comparison_label = comparison.get_pair();

            // ---- Extract Vector of pedigrees
            let pedigree_vec = self.get_pedigree_vec(&comparison_label).loc(loc_msg)?;

            let mut writer = GenericWriter::new(Some(output_files[&comparison_label].clone())).loc(loc_msg)?;
            writer.write_iter(vec![&simulations_header])?;
            writer.write_iter(vec![&pedigree_vec])?;

            if log::log_enabled!(log::Level::Trace) {
                eprintln!("{comparison_label:-^152}");
                eprintln!("{simulations_header}");
                eprintln!("{pedigree_vec}");
            }
        }
        Ok(())
    }

    fn get_most_likely_relationship_zscore(&self, observed_avg_pwd: f64, stats: &[(String, (f64, f64))]) -> (Option<String>, f64, Vec<f64>) {
        let mut min_z_score         : f64            = f64::MAX;
        let mut most_likely_avg_pwd : f64            = 0.0;
        let mut most_likely_rel     : Option<String> = None;

        let mut z_scores = Vec::new();

        for (scenario, (avg_avg_pwd, std_dev)) in stats.iter().rev(){
            let scenario_z_score = (avg_avg_pwd - observed_avg_pwd) / std_dev;
            z_scores.push(scenario_z_score);
            let is_more_likely   = scenario_z_score.abs() < min_z_score.abs();
            if is_more_likely {
                min_z_score         = scenario_z_score;
                most_likely_rel     = Some(scenario.to_owned());
                most_likely_avg_pwd = *avg_avg_pwd;
            }
        }

        (most_likely_rel, most_likely_avg_pwd, z_scores)
    }

    fn get_most_likely_relationship_svm(&self, comparison: &Comparison, pedigree_vec: &PedigreeReps, ordered_rels: &[&String]) -> Result<(Option<String>, Vec<f64>)> {


        let mut svm_builder = LibSvmBuilder::default();
        svm_builder.sims(pedigree_vec).label_order(ordered_rels).scale()?;
        
        let mut predictors = Vec::with_capacity(ordered_rels.len());

        let mut svm_probs  = Vec::new();
        debug!("Estimating most likely relationship through SVM classification");
        for (i, scenario) in ordered_rels[0..ordered_rels.len() - 1].iter().enumerate() {
            let loc_msg = || format!("While attempting to assess likelihood of comparison '{scenario}' for {}", comparison.get_pair());

            predictors.push(ManuallyDrop::new(svm_builder.labels(pedigree_vec, &i).build()?));
            let prediction = predictors[i].predict_probs(comparison.get_avg_pwd()).with_loc(loc_msg)?;
            let true_label = predictors[i].true_label_idx().with_loc(loc_msg)?;
            
            svm_probs.push(prediction[0].1[true_label]);
            if log::log_enabled!(log::Level::Debug){
                println!("  - {i}, P(k>{scenario}): {prediction:?} (index: {true_label})");
            }
        }

        unsafe { predictors.iter_mut().for_each(|svm| ManuallyDrop::drop(svm)) };
        drop(svm_builder);

        // ---- Compute per-class SVM probabilities
        let mut per_class_svm_prob = Vec::with_capacity(svm_probs.len());

        per_class_svm_prob.push(1.0 - svm_probs[0]);                            // k = 1
        for i in 1..svm_probs.len() {
            per_class_svm_prob.push(svm_probs[i-1] - svm_probs[i]);             // 1 < k < K
        }
        per_class_svm_prob.push(svm_probs[svm_probs.len() - 1]);                // k = K

    let max_prob_index = per_class_svm_prob.iter().enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .ok_or(anyhow!("Failed to obtain index of maximum SVM probability"))
        .unwrap();

        let most_likely_rel = ordered_rels.get(max_prob_index).map(|s| (*s).clone());
        Ok((most_likely_rel, per_class_svm_prob))
    }

    /// Recompute a corrected average-PWD of our real samples, and estimate the most likely relationship
    /// using our simulation results
    /// # Arguments
    /// - `comparisons`   : pileup Comparisons of our real samples.
    /// - `output_files`  : target output file where results are written.
    pub fn compute_results(&self, comparisons: &mut PileupComparisons, output_file: &str, assign_method: RelAssignMethod) -> Result<()> {

        let loc_msg = "While attempting to compute corrected summary statistics";
        // ---- Resource acquisition
        let mut simulations_results = Vec::with_capacity(comparisons.len() + 1);
        let mut svm_results         = Vec::with_capacity(comparisons.len() + 1);
        let mut writer = GenericWriter::new(Some(output_file)).loc(loc_msg)?;

        // ---- Print header and write to output_file
        let simulation_header = format!("{: <20} - {: <20} - {: <10} - {: <10} - {: <10} - {: <12} - {: <14} - {: <10} - {: <10}",
            "Pair_name", "Most_Likely_rel", "Corr.Overlap", "Corr.Sum.PWD", "Corr.Avg.PWD", "Corr.CI.95", "Corr.Avg.Phred", "Sim.Avg.PWD", "Min.Z_Score"
        );
 
        simulations_results.push(simulation_header);

        // ---- We might have removed some PWDs during simulations. We need to recompute the variance before printing out
        //      "corrected" 95% Confidence intervals.
        //      -> Run two-pass variance estimation algorithm.
        comparisons.update_variance_unbiased();

        // The use of this indexer helps us ensure SVM columns always keep the same order within the .prob file
        let mut svm_indexer         = IndexMap::with_capacity(comparisons.len());

        // ---- loop across our pileup comparisons and assign a most-likely relationship, using our simulations.
        for comparison in comparisons.iter() {
            let observed_avg_pwd    = comparison.get_avg_pwd();
            let comparison_label    = comparison.get_pair();

            // ---- Gain access to pedigree vector.
            //      @TODO: This error handling is performed multiple times. Stay DRY & wrap this in a public method.
            let pedigree_vec = self.get_pedigree_vec(&comparison_label).loc(loc_msg)?;
    
            // ---- Aggregate the sum of avg. PWD for each relatedness scenario.
            let stats = pedigree_vec.compute_sum_simulated_stats()?;

            let ordered_rels: Vec<&String> = stats.iter().rev().map(|(rel, _)| rel ).collect();

            // ---- Select the scenario having the least amount of Z-score with our observed avg.PWD
            let (mut most_likely_rel, most_likely_avg_pwd, z_scores) = self.get_most_likely_relationship_zscore(observed_avg_pwd, &stats);
            
            // ---- Get svm probabilities if this was assigned.
            if matches!(assign_method, RelAssignMethod::SVM) {
                let (most_likely_rel_svm, svm_probs) = self.get_most_likely_relationship_svm(comparison, pedigree_vec, &ordered_rels)?;
                most_likely_rel = most_likely_rel_svm;

                // ---- Insert label wise svm probabilities for that comparison.
                for i in 0..ordered_rels.len() {
                    *svm_indexer.entry(ordered_rels[i].clone()).or_insert(f64::NAN) = svm_probs[i];
                }

                // ---- Print SVM results to shell if debug mode...
                if log::log_enabled!(log::Level::Debug) {
                    debug!("SVM results for {comparison_label}");
                    println!("  {:<20}{:>12}{:>20}", "Label", "Z-score", "per-class SVM prob");
                    for i in 0..ordered_rels.len() {
                        println!("  {:<20}{:>12.7}{:>20.7}", ordered_rels[i], z_scores[i], svm_probs[i]);
                    }
                    println!("Most likely relationship: {}", &most_likely_rel.as_deref().unwrap_or("None"));
                }
            }
        
            let assigned_rel = most_likely_rel.unwrap_or_else( || {
                warn!("Failed to assign a most likely relationship for {comparison}");
                String::from("None")
            });

            // ---- Get minimum z-score, using the position of ordered_rels and the assigned_rel
            let mut min_z_score = f64::NAN;

            if let Some(min_z_score_index) = ordered_rels.iter().position(|&r| r == &assigned_rel) {
                 min_z_score = z_scores[min_z_score_index]
            }

            
            // ---- Get the "corrected" observed pwd summary statistics. 
            let corrected_sum_pwd = comparison.get_sum_pwd();
            let corrected_overlap = comparison.get_overlap();
            let corrected_ci      = comparison.get_confidence_interval();
            let corrected_phred   = comparison.get_avg_phred();

            // ---- Preformat and log result to console.
            //"Pair_name", "Most_Likely_rel", "Corr.Overlap", "Corr.Sum.PWD", "Corr.Avg.PWD", "Corr.CI.95", "Corr.Avg.Phred", "Sim.Avg.PWD", "Min.Z_Score", 
            let simulation_result = format!(
                "{comparison_label: <20} - \
                 {assigned_rel: <20} - \
                 {corrected_overlap: <12.6} - \
                 {corrected_sum_pwd: <12.6} - \
                 {observed_avg_pwd: <12.6} - \
                 {corrected_ci: <12.6} - \
                 {corrected_phred: <14.6} - \
                 {most_likely_avg_pwd: <11.6} - \
                 {min_z_score: >11.6}"
            );
            simulations_results.push(simulation_result);

            let mut svm_row = format!("{comparison_label:<20} - {observed_avg_pwd:<12.6}");

            for (_, prob) in svm_indexer.iter() {
                svm_row.push_str(&format!(" - {prob:>12.6}"));
            }
            svm_results.push(svm_row);

        }

        // ---- Print simulation results to shell if debug mode.
        if log::log_enabled!(log::Level::Debug) {
            debug!("Simulation results:");
            simulations_results.iter().for_each(|row| println!("{row}"));
        }

        // ---- Write simulation results to file...
        info!("Writing simulation results to output file...");
        writer.write_iter(simulations_results).loc(loc_msg)?;

        // --- Write per-class SVM probability results to output file.
        if matches!(assign_method, RelAssignMethod::SVM) {
            info!("Writing per-class SVM probabilities to output file...");
            let mut svm_results_header = format!("{:<20} - {: <10}", "Pair_name", "Corr.Avg.PWD");
            svm_indexer.keys().for_each(|rel| svm_results_header += &format!(" - {rel: <10}") );

            let svm_output_file   = output_file.strip_suffix("result").unwrap().to_owned() + "probs";
            let mut svm_writer        = GenericWriter::new(Some(svm_output_file)).loc(loc_msg)?;
            svm_writer.write_iter([svm_results_header].into_iter().chain(svm_results)).loc(loc_msg)?;
        }        

        Ok(())
    }
}
