use log::{info, debug, warn};
use located_error::prelude::*;
//pub mod io;
use grups_io::{
    parse::{self, FileKey},
    read::PanelReader,
    read::genotype_reader::{GenotypeReader, VCFReader, FSTReader},
};

mod svm;

pub mod pedigrees;

use pwd_from_stdin::comparisons::Comparisons;


// @TODO! MAIN
//   + [FEATURE] Add ability for multiple contaminating pop
//   + [FEATURE] Add weighted averaged contamination rates ? (see grups_module.pyx:585):
// -------------------------------------------------------------------------------------------------------------------
//   + [  QoL  ] Add multiple candidate file extensions to GeneticMap.
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! CLEANUP + BUGFIX
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! META
// -------------------------------------------------------------------------------------------------------------------
//
pub fn run(
    com_cli           : &parser::Common,
    ped_cli           : &parser::PedigreeSims,
    requested_samples : &[usize],
    comparisons       : &mut Comparisons,
) -> Result<()>
{
    info!("Running 'pedigree-sims' module...");
    // ----------------------------- Sanity checks 
    if ped_cli.contamination_rate.len() < requested_samples.len() {
        warn!("Number of provided contamination rates is lower than that of --samples. \
            Values will be recycled."
        );
    }

    match &ped_cli.seq_error_rate {
        None => {
            // Explicitely warn the user that contamination error rates will be taken from the pileup file
            // if --seq_error_rate was unspecified
            warn!("--seq_error_rate was unspecified. Error probabilities will be sampled directly from the pileup file" );
        },
        Some(error_rates_vec) => {
            if error_rates_vec.len() < requested_samples.len() {
                warn!("Number of provided sequencing error rates is lower than that of --samples. \
                    Values will be recycled."
                );
            }
        }
    }

    // ----------------------------- Prepare output files
    // ---- Add final_results files.
    let mut output_files = parse::get_output_files(
        &mut com_cli.get_file_prefix(None)?, // extract the user requested file prefix
        com_cli.overwrite,                   // Should we allow file overwriting ?
        FileKey::Ext,                        // What key are we using to hash these files ?
        &["".to_string()],                   // Vector of filename suffixes.
        &["result"]                          // Vector of file extensions.
    )?;

    // ---- Add simulations files.
    output_files.extend(
        parse::get_output_files(
            &mut com_cli.get_file_prefix(Some("simulations/"))?,
            com_cli.overwrite,
            FileKey::Suffix,
            &comparisons.get_pairs(),
            &["sims"]
        )?);

    debug!("Output files: {:#?}", output_files);

    // --------------------- Fetch the input panel.
    let mut panel = match ped_cli.panel.as_ref() {
        Some(path) => PanelReader::new(path),
        None       => PanelReader::from_dir(&ped_cli.data_dir),
    }?;

    let input_paths = match ped_cli.mode {
        parser::Mode::Vcf => {
            // --------------------- Get the list of input vcfs.
            info!("Fetching input VCF files in {}", &ped_cli.data_dir.to_str().unwrap_or("None"));
            let input_vcf_paths = VCFReader::fetch_input_files(&ped_cli.data_dir)?;
            panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;
            input_vcf_paths
        },
        parser::Mode::Fst    => FSTReader::<Vec<u8>>::fetch_input_files(&ped_cli.data_dir)?,
        parser::Mode::FstMmap => FSTReader::<memmap2::Mmap>::fetch_input_files(&ped_cli.data_dir)?
    };

    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    info!("Initializing pedigree replicates...");
    let mut pedigrees = pedigrees::Pedigrees::initialize(
        &ped_cli.pedigree_pop,
        comparisons,
        &ped_cli.recomb_dir
    )?;

    // -------------------- Populate all pedigree replicates.
    info!("Populating pedigree replicates...");
    pedigrees.populate(
        comparisons,
        &panel,
        ped_cli.reps,
        &ped_cli.pedigree,
        &ped_cli.contam_pop,
        &ped_cli.contam_num_ind,
    )?;

    // --------------------- Sanity check: if x-chromosome-mode, sex should be defined by the pedigree at this point
    if com_cli.x_chromosome_mode && ! pedigrees.all_sex_assigned() {
        return Err(anyhow!(
            "Using --x-chromosome-mode without defining sexes within the pedigree definition file will produce invalid results. \
            For X-chromosome kinship analysis, please provide the software with a pedigree definition file containing sex assignments for every individuals. \
            Check the documentation for additional information."
        )).loc("While ensuring proper sex assignment of pedigree individuals before X-chromosome analysis.")
    } 

    // --------------------- Randomly assign chromosomal sex of samples if requested
    if com_cli.sex_specific_mode {
        pedigrees.assign_random_sex().loc("While attempting to randomly assign sexes of all pedigrees")?
    }

    // -------------------- Fetch and assign reference sample tags in panel for all founders
    pedigrees.set_founder_tags(&panel).loc("While attempting to randomly assign founder tags of founder individuals in pedigrees")?;  

    // -------------------- Randomly assign and initialize strand-provenance of offpsrings.
    pedigrees.assign_offspring_strands().loc("While attempting to randomly assign strand provenance of offpsrings in pedigrees")?;

    // --------------------- Assign simulation parameters for each pedigree.
    info!("Assigning simulation parameters...");
    pedigrees.set_params(
        comparisons,
        ped_cli.snp_downsampling_rate,
        ped_cli.af_downsampling_rate,
        &ped_cli.seq_error_rate, 
        &ped_cli.contamination_rate
    )?;

    
    // --------------------- Perform pedigree simulations for each pedigree, using all chromosomes.
    match ped_cli.mode {
        parser::Mode::Vcf => {
            info!("Starting VCF pedigree comparisons.");
            for vcf in &input_paths {
                pedigrees.pedigree_simulations_vcf(comparisons,  vcf, ped_cli.maf, ped_cli.decompression_threads)?;
            }
        },
        parser::Mode::Fst => {
            info!("Starting FST pedigree comparisons (in RAM).");
            for fst in &input_paths{
                pedigrees.pedigree_simulations_fst::<Vec<u8>>(comparisons, fst, ped_cli.maf)?;
            }
        },
        parser::Mode::FstMmap => {
            info!("Starting FST pedigree comparisons (Memmap).");
            for fst in &input_paths{
                pedigrees.pedigree_simulations_fst::<memmap2::Mmap>(comparisons, fst, ped_cli.maf)?;
            }
        }
    }

    // --------------------- Print pedigree simulation results.
    pedigrees.write_simulations(comparisons, &output_files)?;


    // --------------------- Compute most likely relationship for each Comparison
    info!("Assigning most likely relationships using {}...", ped_cli.assign_method);
    pedigrees.compute_results(comparisons, &output_files["result"], ped_cli.assign_method)?;

    Ok(())
}