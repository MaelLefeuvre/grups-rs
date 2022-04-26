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

pub fn run(
    //_com_cli          : parser::Common,
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
    let panel =io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path()).unwrap();

    // --------------------- Set ThreadPool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(ped_cli.threads)
        .build()
        .unwrap();


    let template_pedigree= io::pedigree_parser(ped_cli.pedigree.as_path(), &genome).unwrap();
    let mut pedigrees = Vec::new();
    println!("Hi.");
    for i in 0..500 {
        let mut new_pedigree = template_pedigree.clone();
        new_pedigree.set_tags(&panel, &"EUR".to_string());
        pedigrees.push(new_pedigree);
    }
    println!("Done");



    for comparison in comparisons.get() {
        pool.scope(|scope| {
            for _ in 0..ped_cli.reps {
                scope.spawn(|_| {
                    let thread_idx = pool.current_thread_index().unwrap();
                    debug!("[Thread {thread_idx}]: spawned!", );

                    // --------------------- Parse input pedigree File
                    info!("[Thread {thread_idx}]: Parsing template pedigree.");

                    let mut pedigree= io::pedigree_parser(ped_cli.pedigree.as_path(), &genome).unwrap();
                    pedigree.set_tags(&panel, &"EUR".to_string());
                    for founder in pedigree.founders_mut() {
                        info!("[Thread {thread_idx}]: Founder   : {}", founder);
                    }
                    for offspring in pedigree.offsprings_mut() {
                        info!("[Thread {thread_idx}]: Offspring : {}", offspring);
                    }
                    pedigree.populate_founders_vcf(&input_vcf_paths, &comparison.positions).unwrap();
                    pedigree.reproduce(&genetic_map);
                    pedigree.compare_genomes();
                });
            }
        });
    }
    Ok(())
}