pub mod pileup;
pub mod comparison;
pub mod io;

//extern crate logger;
//use parser;

use comparison::*;
use genome::{
    SNPCoord,
    Genome,
};
use io::*;

use std::fs;
use std::error::Error;
use std::collections::HashSet;
use std::io::{BufReader, BufRead};

use log::{error, warn, info, debug};
use std::process;


// Main function. 
pub fn run<'a>(
    com_cli           : &'a parser::Common,
    pwd_cli           : &'a parser::PwdFromStdin,
    requested_samples : &'a [usize],
    genome            : &'a Genome,
) -> Result<(Comparisons, HashSet<SNPCoord>), Box<dyn Error>>{

    // ----------------------------- Sanity checks.
    pwd_cli.check_depth()?; // Ensure min_depths are > 2 when allowing self-comparisons
    com_cli.check_input()?; // Ensure the user has either requested stdin or --pileup

    if pwd_cli.min_depth.len() < requested_samples.len() {
        warn!("--min-depth length is less than that of --samples. Values of min-depth will wrap around.")
    }

    // ----------------------------- Parse Comparisons
    info!("Parsing Requested comparisons...");
    let mut comparisons = Comparisons::parse(requested_samples, &pwd_cli.min_depth, &com_cli.sample_names, pwd_cli.self_comparison, genome, pwd_cli.blocksize);

    // ----------------------------- Prepare output files
    // ---- Add pwd files.
    let mut output_files = io::get_output_files(
        &mut com_cli.get_file_prefix(None).unwrap(),    // extract the user requested file prefix
        com_cli.overwrite,                                  // Should we allow file overwriting ?
        FileKey::Ext,         // What key are we using to hash these files ?
        &["".to_string()],// Vector of filename suffixes.
        &["pwd"]          // Vector of file extensions.
    )?;

    // ---- Add blocks files.
    output_files.extend(
        io::get_output_files(
            &mut com_cli.get_file_prefix(Some("blocks/")).unwrap(),
            com_cli.overwrite,
            FileKey::Suffix,
            &comparisons.get_pairs(),
            &["blk"]
        )?.into_iter());

    debug!("Output files: {:#?}", output_files);


    // ----------------------------- Parse target_positions
    info!("Parsing target_positions...");
    let target_positions = match &com_cli.targets {
        None => HashSet::new(),
        Some(filename) => SNPReader::new(filename)?.hash_target_positions()?
    };

    let target_required: bool = ! target_positions.is_empty();

    // ----------------------------- Parse requested Chromosomes
    let valid_chromosomes : Vec<u8> = match com_cli.chr.clone() {
        None         => genome.keys().copied().collect(),
        Some(vector) => parser::parse_user_ranges(&vector, "chr")?
    };
    info!("Valid chromosomes: {:?}", valid_chromosomes);

    // ---------------------------- Choose between file handle or standard input
    info!("Opening pileup...");   
    let pileup_reader: Box<dyn BufRead> = match &com_cli.pileup {
        None => Box::new(BufReader::new(std::io::stdin())),
        Some(filename) => Box::new(BufReader::new(fs::File::open(filename).unwrap()))
    };

    // ---------------------------- Read Pileup
    info!("Parsing pileup...");   
    for entry in pileup_reader.lines() {
        // ----------------------- Parse line.
        let mut line: pileup::Line = match pileup::Line::new(entry.as_ref().unwrap(), pwd_cli.ignore_dels){
            Ok(line) => line,
            Err(e) => {error!("Error: {}", e); process::exit(1);},
        };

        // ----------------------- Check if line should be skipped.
        if ! valid_chromosomes.contains(&line.coordinate.chromosome) {
            continue  // Skip line if this is not a valid chromosome. 
        }

        if target_required && !target_positions.contains(&line.coordinate) {  
            continue  // Skip line if we're targeting snps + the current line is not found. 
        }

        // ------------------------ Apply quality filtering on all individuals.
        line.filter_base_quality(&com_cli.min_qual);

        // ------------------------ Apply target filtration if requested.
        if pwd_cli.known_variants {
            let current_coord = match target_positions.get(&line.coordinate) {
                Some(coordinate) => coordinate,
                None => panic!("Cannot filter known variants when REF/ALT allele are unknown! Please use a different file format."),
            };
            line.filter_known_variants(current_coord);
        }

        // ----------------------- Compute PWD (or simply print the line if there's an existing overlap)        
        for comparison in comparisons.get_mut() {
            if comparison.satisfiable_depth(&line.individuals) {
                if ! pwd_cli.filter_sites {
                    comparison.compare(&line)?;
                } else {
                    println!("{}", entry.as_ref().unwrap());
                }
            }
        }
    }
    

    // ----------------------------- Print results
    let mut pwd_writer = Writer::new(Some(output_files["pwd"].clone()))?;
    if ! pwd_cli.filter_sites {
        info!("Printing results...");
        let header = format!("{: <20} - Overlap - Sum PWD - Avg. Pwd - Avg. Phred", "Name");
        println!("{}", header);
        pwd_writer.write_iter(&vec![header])?;       // Print PWD results to file.
        pwd_writer.write_iter(comparisons.get())?;   // 

        println!("{}", comparisons);                 // Print PWD results to console

        if pwd_cli.print_blocks {
            for comparison in comparisons.get() {
                let pair = comparison.get_pair();
                let mut block_writer = Writer::new(Some(output_files[&pair].clone()))?;
                block_writer.write_iter(vec![&comparison.blocks])?;
            }
        }
    }
    Ok((comparisons, target_positions))
}
