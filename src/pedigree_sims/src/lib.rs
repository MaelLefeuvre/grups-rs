use std::{
    error::Error, collections::HashMap, path::PathBuf,
};
use log::{info};
//use rayon;

use genome::{
    Genome,
    GeneticMap,
};
use pwd_from_stdin::comparison::Comparisons;
pub mod io;
pub mod pedigree;

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
    // --------------------- Parse input recombination maps.
    info!("Parsing genetic maps in {}", &ped_cli.recomb_dir.to_str().unwrap_or("None"));
    let genetic_map = GeneticMap::default().from_dir(&ped_cli.recomb_dir)?;

    // --------------------- Fetch the input panel.
    let panel = match ped_cli.panel.clone() {
        Some(path) => path,
        None => io::vcf::fetch_input_panel(&ped_cli.data_dir)?,
    };

    // --------------------- Parse Input Samples Panel
    let mut panel = io::vcf::reader::VCFPanelReader::new(panel.as_path())?;

    let vcf_requested = false;

    let input_paths = if vcf_requested {
        // --------------------- Get the list of input vcfs.
        info!("Fetching input VCF files in {}", &ped_cli.data_dir.to_str().unwrap_or("None"));
        let input_vcf_paths = io::vcf::get_input_vcfs(&ped_cli.data_dir)?;
        panel.assign_vcf_indexes(input_vcf_paths[0].as_path())?;
        input_vcf_paths
    }
    else {
        let fst_dir = PathBuf::from("./tests/test-data/fst/");
        let input_fst_paths = io::fst::get_input_fst(&fst_dir)?;
        input_fst_paths
    };


    // --------------------- Get random contaminating individuals. [PROTOTYPE]
    let mut contam_ind_ids = Vec::new();
    for _ in 0..1 {
        let contam_sample_tag = panel.random_sample(&ped_cli.contam_pop[0]).unwrap();
        contam_ind_ids.push(contam_sample_tag);
    }

    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    let template_pedigree= io::pedigree::pedigree_parser(ped_cli.pedigree.as_path(), &genome).unwrap();
    let mut pedigrees = HashMap::new();
    for comparison in comparisons.get() {
        let comparison_label = comparison.get_pair();
        pedigrees.insert(comparison_label.to_owned(), Vec::new());
        for _ in 0..ped_cli.reps {
            let mut new_pedigree = template_pedigree.clone();
            new_pedigree.set_tags(&panel, &ped_cli.pedigree_pop);
            new_pedigree.assign_offspring_strands()?;
            pedigrees.get_mut(&comparison_label).unwrap().push(new_pedigree);
        }
    }

    // --------------------- Set ThreadPool
    //let pool = rayon::ThreadPoolBuilder::new()
    //    .num_threads(ped_cli.threads)
    //    .build()
    //    .unwrap();

    // --------------------- Perform pedigree simulations for each pedigree, using all chromosomes.
    if vcf_requested {
        info!("Starting VCF pedigree comparisons.");
        for vcf in input_paths.iter() {
            let simulations = pedigree::pedigree_simulations(
                &mut pedigrees, 
                vcf,
                comparisons,
                &ped_cli.pedigree_pop,
                ped_cli.contamination_rate[0][0] as f64 / 100.0,
                &contam_ind_ids,
                ped_cli.pmd_rate[0][0] as f64 / 100.0,
                ped_cli.af_downsampling_rate,
                ped_cli.snp_downsampling_rate,
                &genetic_map,
                ped_cli.maf,
                ped_cli.decompression_threads
            );
            simulations?;
        }
    } else {
           // -------------- [PROTOTYPE]
        info!("Starting FST pedigree comparisons.");
        for fst in input_paths.iter(){
            let simulations = pedigree::pedigree_simulations_fst(
                &mut pedigrees, 
                fst,
                comparisons,
                &ped_cli.pedigree_pop,
                ped_cli.contamination_rate[0][0] as f64 / 100.0,
                &contam_ind_ids,
                ped_cli.pmd_rate[0][0] as f64 / 100.0,
                ped_cli.af_downsampling_rate,
                ped_cli.snp_downsampling_rate,
                &genetic_map,
                ped_cli.maf,
                ped_cli.decompression_threads
            );
            simulations?;
        }
        // -------------- [END PROTOTYPE]
    }

    // --------------------- Print pedigree simulation results.
    for req_comparison in comparisons.get() {
        let comparison_label = req_comparison.get_pair();
        println!("--------------------- {comparison_label}");
        let pedigree_vec = pedigrees.get_mut(&comparison_label).unwrap();
        for (i, pedigree) in pedigree_vec.iter().enumerate() {
            pedigree.print_results(i);
        }
    }

    // --------------------- Compute most likely relationship for each Comparison
    println!("----------------------------------------------------");
    for req_comparison in comparisons.get() {
        let comparison_label = req_comparison.get_pair();
        let pedigree_vec = pedigrees.get_mut(&comparison_label).unwrap();
        pedigree::compute_results(&pedigree_vec, req_comparison);
    }


    println!("Done!");

    Ok(())
}