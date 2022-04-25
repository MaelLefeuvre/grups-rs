use std::{
    collections::HashSet,
    error::Error,
};
use log::{info};
use pwd_from_stdin::genome::{SNPCoord, Genome, self};
use pwd_from_stdin::comparison::Comparisons;
//use parser;
pub mod io;
pub mod pedigree;

//use rust_htslib::tpool::ThreadPool;


pub fn run<'a>(
    _com_cli          : &'a parser::Common,
    ped_cli           : &'a parser::PedigreeSims,
    _requested_samples: &'a [usize],
    genome            : &'a Genome,
    comparisons       : &'a Option<Comparisons>,
    _target_positions : &'a Option<HashSet<SNPCoord>>
) -> Result<(), Box<dyn Error>>
{

    let tpool = rust_htslib::tpool::ThreadPool::new(4).unwrap();

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
    let panel = io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path(), &tpool).unwrap();
    //println!("{:?}", panel.samples[&"EUR".to_string()]);

    // --------------------- Parse input pedigree File
    info!("Parsing template pedigree.");
    let mut pedigree= io::pedigree_parser(ped_cli.pedigree.as_path(), genome).unwrap();
    for founder in pedigree.founders_mut() {
        info!("Founder   : {}", founder);
    }

    for offspring in pedigree.offsprings_mut() {
        info!("Offspring : {}", offspring);
    }

    // TODO: Assign contaminant to pedigree template

    for comparison in comparisons.as_ref().unwrap().get() {
        // Multithread starts here.
        for i in 0..4 {
            std::thread::spawn(move || {
                println!("hi number {} from the spawned thread!", i);
            });
        }
        pedigree.set_tags(&panel, &"EUR".to_string());
        pedigree.populate_founders_vcf(&input_vcf_paths, &comparison.positions).unwrap();
        pedigree.reproduce(&genetic_map);
        pedigree.compare_genomes();
        //println!("{:?}", pedigree);
    }
    Ok(())
}