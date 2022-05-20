use std::error::Error;
use log::{info, debug, warn};

pub mod io;
pub mod pedigrees;

use pwd_from_stdin::comparisons::Comparisons;


// @TODO! MAIN
//   +        [CRUCIAL] Add ability for multiple contaminating pop
//   +        [CRUCIAL] Add ability to compute error rate from pileup file (see grups_module.pyx:610). 
//   +                 - comparisons positions should keep a record of phred scales.
//   +        [CRUCIAL] Add weighted averaged contamination rates (see grups_module.pyx:585):
// -------------------------------------------------------------------------------------------------------------------
//   +        [  QoL  ] Finish implementing CLI Args (de)serialization.
//
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! CLEANUP + BUGFIX
//   + [CRUCIAL] Remove dead arguments from CLI parser
//   + [  BUG  ] GeneticMap does not detect if recomb-dir is empty.
//   + [  BUG  ] GeneticMap panics if file != ".txt"
//   + [  BUG  ] Grups Fst bugs out when Samples Panel does not perfectly match the VCF ordering.
//   + [  BUG  ] VcfPanelReader::copy_from_source not working when using FST pop-subset
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! META
//   + [CRUCIAL] Refactor pwd_from_stdin::io and pedigree_sims::io into a self-contained library.
//   + [CRUCIAL] Document pedigree_sims::* libraries.
//
pub fn run(
    _com_cli          : parser::Common,
    ped_cli           : parser::PedigreeSims,
    comparisons       : &Comparisons,
) -> Result<pedigrees::Pedigrees, Box<dyn Error>>
{

    // ----------------------------- Sanity checks 
    if ped_cli.contamination_rate.len() < comparisons.len() {
        warn!("Number of contamination rates is lower than the number of comparisons. Values of --contamination-rate will wrap around.")
    }

    if ped_cli.pmd_rate.len() < comparisons.len() {
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
    let mut pedigrees = pedigrees::Pedigrees::initialize(ped_cli.pedigree_pop, comparisons, &ped_cli.recomb_dir)?;
    pedigrees.populate(&panel, ped_cli.reps, &ped_cli.pedigree, &ped_cli.contam_pop, &ped_cli.contam_num_ind)?;

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

    Ok(pedigrees)
}