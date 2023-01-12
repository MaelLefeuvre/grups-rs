use log::{info, debug, warn};
use located_error::prelude::*;
//pub mod io;
use grups_io::{
    parse::{self, FileKey},
    read::PanelReader,
    read::genotype_reader::{GenotypeReader, VCFReader, FSTReader},
};

pub mod pedigrees;

use pwd_from_stdin::comparisons::Comparisons;


// @TODO! MAIN
//   +        [CRUCIAL] Add ability for multiple contaminating pop
//   +        [CRUCIAL] Add weighted averaged contamination rates (see grups_module.pyx:585):
//   +        [CRUCIAL] Add Missing SNP filtration mechanism. (remove &Pwd if Pwd.coordinate is not in vcf/fst file.)
//                      - add "typed" counter in Pwd field ? -> remove if typed == 0 after all files have been used.
// -------------------------------------------------------------------------------------------------------------------
//            [  QoL  ] Add multiple candidate file extensions to GeneticMap.
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! CLEANUP + BUGFIX
//   + [  BUG  ] Grups Fst bugs out when Samples Panel does not perfectly match the VCF ordering.
//   + [  BUG  ] VcfPanelReader::copy_from_source not working when using FST pop-subset
//   + [  BUG  ] vcfreader.rs:193:39 && fstreader.rs:144:14 panics if contaminating pop is missing from vcf/fst index.
// -------------------------------------------------------------------------------------------------------------------
// @ TODO! META
//   +          [CRUCIAL] Refactor pwd_from_stdin::io and pedigree_sims::io into a self-contained library.
//   + [ DONE  ][CRUCIAL] Document pedigree_sims::* libraries.
//
pub fn run(
    com_cli           : parser::Common,
    ped_cli           : parser::PedigreeSims,
    comparisons       : &mut Comparisons,
) -> Result<()>
{

    // ----------------------------- Sanity checks 
    if ped_cli.contamination_rate.len() < comparisons.len() {
        warn!("Number of contamination rates is lower than the number of comparisons. \
            Values of --contamination-rate will wrap around."
        );
    }

    match &ped_cli.seq_error_rate {
        None => {
            // Explicitely warn the user that contamination error rates will be taken from the pileup file
            // if --seq_error_rate was unspecified
            warn!("--seq_error_rate was unspecified. Error probabilities will be sampled directly from the pileup file" );
        },
        Some(error_rates_vec) => {
            if error_rates_vec.len() < comparisons.len() {
                warn!("Number of sequencing error rates is lower than the number of comparisons. \
                    Values of --seq-error-rate will wrap around."
                );
            }
        }
    }

    // ----------------------------- Prepare output files
    // ---- Add final_results files.
    let mut output_files = parse::get_output_files(
        &mut com_cli.get_file_prefix(None)?, // extract the user requested file prefix
        com_cli.overwrite,                       // Should we allow file overwriting ?
        FileKey::Ext,                   // What key are we using to hash these files ?
        &["".to_string()],   // Vector of filename suffixes.
        &["result"]          // Vector of file extensions.
    )?;

    // ---- Add blocks files.
    output_files.extend(
        parse::get_output_files(
            &mut com_cli.get_file_prefix(Some("simulations/"))?,
            com_cli.overwrite,
            FileKey::Suffix,
            &comparisons.get_pairs(),
            &["sims"]
        )?.into_iter());

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
        parser::Mode::Fst => {
            FSTReader::fetch_input_files(&ped_cli.data_dir)?
        },
    };

    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    info!("Initializing pedigree replicates...");
    let mut pedigrees = pedigrees::Pedigrees::initialize(
        ped_cli.pedigree_pop,
        comparisons,
        &ped_cli.recomb_dir
    )?;

    info!("Populating pedigree replicates...");
    pedigrees.populate(
        comparisons,
        &panel,
        ped_cli.reps,
        &ped_cli.pedigree,
        &ped_cli.contam_pop,
        &ped_cli.contam_num_ind
    )?;

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
                pedigrees.pedigree_simulations_vcf(
                    comparisons, 
                    vcf,
                    ped_cli.maf,
                    ped_cli.decompression_threads
                )?;
            }
        },
        parser::Mode::Fst => {
            info!("Starting FST pedigree comparisons.");
            for fst in &input_paths{
                pedigrees.pedigree_simulations_fst(
                    comparisons,
                    fst,
                    ped_cli.maf
                )?;
            }
        },
    }

    //pedigrees.filter(comparisons);

    // --------------------- Print pedigree simulation results.
    pedigrees.write_simulations(comparisons, &output_files)?;


    // --------------------- Compute most likely relationship for each Comparison
    println!("--------------------------------------------------------------");
    pedigrees.compute_results(comparisons, &output_files["result"])?;

    Ok(())
}