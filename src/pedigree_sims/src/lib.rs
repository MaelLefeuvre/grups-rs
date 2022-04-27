use std::{
    error::Error,
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
    //_requested_samples: &'a [usize],
    genome            : Genome,
    comparisons       : &Comparisons,
    //_target_positions : Option<HashSet<SNPCoord>>
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


    let template_pedigree= io::pedigree_parser(ped_cli.pedigree.as_path(), &genome).unwrap();
    let mut pedigrees = Vec::new();
    for _ in 0..ped_cli.reps {
        let mut new_pedigree = template_pedigree.clone();
        new_pedigree.set_tags(&panel, &ped_cli.pedigree_pop);
        new_pedigree.assign_offspring_strands()?;
        pedigrees.push(new_pedigree);
    }

    // Perform pedigree simulations for each pedigree, using all chromosomes.
    info!("Starting pedigree comparisons.");
    for comparison in comparisons.get() {
        comparison.get_pair();
        for vcf in input_vcf_paths.iter() {
            let simulations = pedigree::pedigree_simulations(&mut pedigrees, vcf, &comparison.positions, &ped_cli.pedigree_pop, &genetic_map, ped_cli.threads);
            simulations.await?;
        }
        // Print pedigree simulation results.
        for req_comparison in comparisons.get() {
            //println!("----------------- {}", req_comparison.get_pair());
            for (i, pedigree) in pedigrees.iter().enumerate() {
                pedigree.print_results(i);
            }
        }
        pedigree::compute_results(&pedigrees, comparison);
    }
    println!("Done!");

    Ok(())
}