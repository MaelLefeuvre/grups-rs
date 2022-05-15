use std::{
    error::Error,
};
use log::{info, debug, warn};
//use rayon;

use genome::{
    Genome
};

pub mod io;
pub mod pedigree;
pub mod pedparam;
pub mod contaminant;
use pwd_from_stdin::comparison::Comparisons;


// @TODO 
//   + [DONE] [CRUCIAL] Print simulations results into files.
//   + [DONE] [CRUCIAL] Pop selected contaminating individuals out of the pedigree founders candidates.
//   +        [CRUCIAL] Add ability to compute error rate from pileup file (see grups_module.pyx:610). 
//   +                 - comparisons positions should keep a record of phred scales.
//   + [DONE] [CRUCIAL] Add weighted averaged contamination rates (see grups_module.pyx:585):
//   + [DONE] [CRUCIAL] Find a way to restrict the number of tokio worker-threads!
// -------------------------------------------------------------------------------------------------------------------
//   + [DONE] [FEATURE] Implement ability to add multiple --contam-pop
//   + [DONE] [FEATURE] ability to add multiple --contam-num-inds
//   + [DONE] [FEATURE] ability to add different contaminating individuals per pair. 
//   + [DONE] [FEATURE] Add nSNP downsampling (ds_rate_num_snps)
//   + [DONE] [FEATURE] Add range ability for contamination-rate, pmd-rate, depth-rate
//   +        [FEATURE] Add built-in vcf filtration module.
// -------------------------------------------------------------------------------------------------------------------
//   +        [ SPEED ] Implement per-vcf parallelization (?)
//                      + Will most probably require refactoring from Rc<RefCell<Individual>> to Arc<Mutex<Individual>>.
//                        + Might prove unproductive, as this could greatly impact memory consumption.
//                        + Still needs some testing, as VCF I/O remains the biggest performance bottleneck.
//   + [DONE] [ SPEED ] Benchmark Finite-State Transducer indexation mode.
// -------------------------------------------------------------------------------------------------------------------
//   + [DONE] [  QoL  ] If multiple parallelization options, separate into --threads and --decompression-threads for more 
//                      flexibility.
//   +        [  QoL  ] Finish implementing CLI Args (de)serialization.
//   + [DONE] [  QoL  ] Add .vcf and GZip compressed support for vcf. (async I/O support for vcf, please...)
// -------------------------------------------------------------------------------------------------------------------
//
// @ TODO! META
//   + [CRUCIAL] Refactor pwd_from_stdin::io and pedigree_sims::io into a self-contained library.
//   + [CRUCIAL] More unit-testing, please!!! 
//   + [CRUCIAL] Document pedigree_sims::* libraries.
//   + [CRUCIAL] Split up pedigree_sims::pedigree::pedigree_simulations() into a managable Object or function, please.
//   + [  QoL  ] Clean up or refactor dead code
//
pub fn run(
    _com_cli          : parser::Common,
    ped_cli           : parser::PedigreeSims,
    genome            : Genome,
    comparisons       : &Comparisons,
) -> Result<(), Box<dyn Error>>
{

    // ----------------------------- Sanity checks 
    if ped_cli.contamination_rate.len() < comparisons.get().len() {
        warn!("Number of contamination rates is lower than the number of comparisons. Values of --contamination-rate will wrap around.")
    }

    if ped_cli.pmd_rate.len() < comparisons.get().len() {
        warn!("Number of sequencing error rates is lower than the number of comparisons. Values of --contamination-rate will wrap around.")
    }

    // ----------------------------- Prepare output files
    // ---- Add final_results files.
    let mut output_files = pwd_from_stdin::io::get_output_files(
        &mut _com_cli.get_file_prefix(None).unwrap(),    // extract the user requested file prefix
        _com_cli.overwrite,                                  // Should we allow file overwriting ?
        pwd_from_stdin::io::FileKey::Ext,         // What key are we using to hash these files ?
        &["".to_string()],   // Vector of filename suffixes.
        &["result"]          // Vector of file extensions.
    )?;

    // ---- Add blocks files.
    output_files.extend(
        pwd_from_stdin::io::get_output_files(
            &mut _com_cli.get_file_prefix(Some("simulations/")).unwrap(),
            _com_cli.overwrite,
            pwd_from_stdin::io::FileKey::Suffix,
            &comparisons.get_pairs(),
            &["sims"]
        )?.into_iter());

    debug!("Output files: {:#?}", output_files);

    // --------------------- Fetch the input panel.
    let panel = match ped_cli.panel.clone() {
        Some(path) => path,
        None => io::vcf::fetch_input_panel(&ped_cli.data_dir)?,
    };

    // --------------------- Parse Input Samples Panel
    let mut panel = io::vcf::reader::VCFPanelReader::new(panel.as_path())?;

    let input_paths = match ped_cli.mode {
        parser::Mode::Vcf => {
            // --------------------- Get the list of input vcfs.
            info!("Fetching input VCF files in {}", &ped_cli.data_dir.to_str().unwrap_or("None"));
            let input_vcf_paths = io::vcf::get_input_vcfs(&ped_cli.data_dir)?;
            panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;
            input_vcf_paths
        },
        parser::Mode::Fst => {
            io::fst::get_input_fst(&ped_cli.data_dir)?
        },
    };

    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    let mut pedigrees = pedigree::Pedigrees::initialize(ped_cli.pedigree_pop, comparisons, &ped_cli.recomb_dir)?;
    pedigrees.populate(&panel, ped_cli.reps, &genome, &ped_cli.pedigree, &ped_cli.contam_pop, &ped_cli.contam_num_ind)?;

    // --------------------- Assign simulation parameters for each pedigree.
    pedigrees.set_params(ped_cli.snp_downsampling_rate, ped_cli.af_downsampling_rate, &ped_cli.pmd_rate, &ped_cli.contamination_rate)?;

    // --------------------- Perform pedigree simulations for each pedigree, using all chromosomes.
    match ped_cli.mode {
        parser::Mode::Vcf => {
            info!("Starting VCF pedigree comparisons.");
            for vcf in input_paths.iter() {
                pedigrees.pedigree_simulations_vcf(
                    vcf,
                    ped_cli.maf,
                    ped_cli.decompression_threads
                )?;
            }
        },
        parser::Mode::Fst => {
            // -------------- [PROTOTYPE]
            info!("Starting FST pedigree comparisons.");
            for fst in input_paths.iter(){
                pedigrees.pedigree_simulations_fst(
                    fst,
                    ped_cli.maf
                )?;
            }
        },
        // -------------- [END PROTOTYPE]
    }

    // --------------------- Print pedigree simulation results.
    pedigrees.write_simulations(&output_files)?;


    // --------------------- Compute most likely relationship for each Comparison
    println!("--------------------------------------------------------------");
    pedigrees.compute_results(&output_files["result"])?;

    Ok(())
}