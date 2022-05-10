use std::{
    error::Error,
};
use log::{info, debug};
//use rayon;

use genome::{
    Genome
};

pub mod io;
pub mod pedigree;

use pwd_from_stdin::comparison::Comparisons;


/// @TODO! 
///   - [CRUCIAL] Print simulations results into files.
///   - [CRUCIAL] Pop selected contaminating individuals out of the pedigree founders candidates.
///   - [CRUCIAL] Add ability to compute error rate from pileup file (see grups_module.pyx:610). 
///               - comparisons positions should keep a record of phred scales.
///   - [CRUCIAL] Add weighted averaged contamination rates (see grups_module.pyx:585):
///   - [CRUCIAL] Find a way to restrict the number of tokio worker-threads!
/// -------------------------------------------------------------------------------------------------------------------
///   - [FEATURE] Implement ability to add multiple --contam-pop
///   - [FEATURE] ability to add multiple --contam-num-inds
///   - [FEATURE] ability to add different contaminating individuals per pair. 
///   - [FEATURE] Add nSNP downsampling (ds_rate_num_snps)
///   - [FEATURE] Add range ability for contamination-rate, pmd-rate, depth-rate
///               + these should be expressed as floats --> itertools::linspace could be a candidate package
///                 to work-out ranges.
///               + param_num_rep might not be the most suitable parameter choice ?
///   - [FEATURE] Add built-in vcf filtration module.
/// -------------------------------------------------------------------------------------------------------------------
///   - [ SPEED ] Implement per-vcf parallelization (?)
///               - Will most probably require refactoring from Rc<RefCell<Individual>> to Arc<Mutex<Individual>>.
///                 => Might prove unproductive, as this could greatly impact memory consumption.
///                 => Still needs some testing, as VCF I/O remains the biggest performance bottleneck.
///   - [ SPEED ] Benchmark Finite-State Transducer indexation mode.
/// -------------------------------------------------------------------------------------------------------------------
///   - [  QoL  ] If multiple parallelization options, separate into --threads and --decompression-threads for more 
///               flexibility.
///   - [  QoL  ] Finish implementing CLI Args (de)serialization.
///   - [  QoL  ] Add .vcf and GZip compressed support for vcf. (async I/O support for vcf, please...)
/// -------------------------------------------------------------------------------------------------------------------
///
/// @ TODO! META
///   - [CRUCIAL] Refactor pwd_from_stdin::io and pedigree_sims::io into a self-contained library.
///   - [CRUCIAL] More unit-testing, please!!! 
///   - [CRUCIAL] Document pedigree_sims::* libraries.
///   - [CRUCIAL] Split up pedigree_sims::pedigree::pedigree_simulations() into a managable Object or function, please.
///   - [  QoL  ] Clean up or refactor dead code
/// 
pub fn run(
    _com_cli          : parser::Common,
    ped_cli           : parser::PedigreeSims,
    genome            : Genome,
    comparisons       : &Comparisons,
) -> Result<(), Box<dyn Error>>
{
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


    // --------------------- Get random contaminating individuals. [PROTOTYPE]
    let mut contam_ind_ids = Vec::new();
    for _ in 0..1 {
        let contam_sample_tag = panel.random_sample(&ped_cli.contam_pop[0]).unwrap();
        contam_ind_ids.push(contam_sample_tag);
    }

    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    let mut pedigrees = pedigree::Pedigrees::new(&ped_cli.pedigree, ped_cli.reps, ped_cli.pedigree_pop, comparisons, &panel, &genome, &ped_cli.recomb_dir)?;

    // --------------------- Set ThreadPool
    //let pool = rayon::ThreadPoolBuilder::new()
    //    .num_threads(ped_cli.threads)
    //    .build()
    //    .unwrap();

    // --------------------- Perform pedigree simulations for each pedigree, using all chromosomes.
    let contam_rate = ped_cli.contamination_rate[0][0] as f64 / 100.0; 
    let seq_error_rate = ped_cli.pmd_rate[0][0] as f64 / 100.0;

    match ped_cli.mode {
        parser::Mode::Vcf => {
            info!("Starting VCF pedigree comparisons.");
            for vcf in input_paths.iter() {
                pedigrees.pedigree_simulations_vcf(
                    vcf,
                    contam_rate,
                    &contam_ind_ids,
                    seq_error_rate,
                    ped_cli.af_downsampling_rate,
                    ped_cli.snp_downsampling_rate,
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
                    contam_rate,
                    &contam_ind_ids,
                    seq_error_rate,
                    ped_cli.af_downsampling_rate,
                    ped_cli.snp_downsampling_rate,
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