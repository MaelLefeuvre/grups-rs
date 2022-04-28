use std::{
    error::Error, collections::HashMap,
};
use log::{info, debug};
//use rayon;

use pwd_from_stdin::{
    genome::{Genome, self},
    comparison::Comparisons,
};

pub mod io;
pub mod pedigree;

use tokio;

#[tokio::main]
pub async fn run(
    _com_cli          : parser::Common,
    ped_cli           : parser::PedigreeSims,
    genome            : Genome,
    comparisons       : &Comparisons,
) -> Result<(), Box<dyn Error>>
{

    // --------------------- Get the list of input vcfs.
    info!("Fetching input VCF files in {}", &ped_cli.data_dir.to_str().unwrap_or("None"));
    let input_vcf_paths = io::get_input_vcfs(&ped_cli.data_dir)?;

    // --------------------- Fetch the input panel.
    let panel = match ped_cli.panel.clone() {
        Some(path) => path,
        None => io::fetch_input_panel(&ped_cli.data_dir)?,
    };

    // --------------------- Parse input recombination maps.
    info!("Parsing genetic maps in {}", &ped_cli.recomb_dir.to_str().unwrap_or("None"));
    let genetic_map = genome::GeneticMap::default().from_dir(&ped_cli.recomb_dir)?;

    // --------------------- Parse Input Samples Panel
    let panel =io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path()).await?;

    // --------------------- Set ThreadPool
    //let pool = rayon::ThreadPoolBuilder::new()
    //    .num_threads(ped_cli.threads)
    //    .build()
    //    .unwrap();


    // --------------------- Generate empty pedigrees for each Comparison & each requested replicate.
    let template_pedigree= io::pedigree_parser(ped_cli.pedigree.as_path(), &genome).unwrap();
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

    // --------------------- Perform pedigree simulations for each pedigree, using all chromosomes.
    info!("Starting pedigree comparisons.");
    for vcf in input_vcf_paths.iter() {
        let simulations = pedigree::pedigree_simulations(
            &mut pedigrees, 
            vcf,
            comparisons,
            &ped_cli.pedigree_pop,
            &genetic_map,
            ped_cli.maf,
            ped_cli.threads
        );
        simulations.await?;
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