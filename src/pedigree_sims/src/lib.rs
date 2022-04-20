use std::{
    collections::HashSet,
    error::Error,
};
use log::{info};
use pwd_from_stdin::genome::{Chromosome, SNPCoord};
use pwd_from_stdin::comparison::Comparisons;
//use parser;
pub mod io;
pub mod pedigree;

pub fn run<'a>(
    _com_cli: &'a parser::Common,
    ped_cli: &'a parser::PedigreeSims,
    _requested_samples: &'a [usize],
    _genome: &'a [Chromosome],
    _comparisons: &'a Option<Comparisons>,
    _target_positions: &'a Option<HashSet<SNPCoord>>
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
    let _recombination_map = io::GenMapReader::new(&ped_cli.recomb_dir)?;
    //println!("recombination_map {:#?}", recombination_map.get(&11, &2190951));

    // --------------------- Parse Input Samples Panel
    let panel = io::VCFPanelReader::new(panel.as_path(), input_vcf_paths[0].as_path()).unwrap();
    //println!("{:?}", panel.samples[&"EUR".to_string()]);

    // --------------------- Parse input pedigree File
    info!("Parsing template pedigree.");
    let mut pedigree= io::pedigree_parser(ped_cli.pedigree.as_path()).unwrap();
    for founder in pedigree.founders_mut() {
        info!("Founder   : {}", founder);
    }

    for offspring in pedigree.offsprings_mut() {
        info!("Offspring : {}", offspring);
    }

    // TODO: Assign contaminant to pedigree template

    // Multithread starts here.
    for mut founder in pedigree.founders_mut() {
            let sample = panel.random_sample(&"EUR".to_string()).unwrap();
            founder.assign_genome(sample, &input_vcf_paths)?;
    }
    
    Ok(())
}