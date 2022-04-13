extern crate pwd_from_stdin;

use crate::pwd_from_stdin::pileup;
use crate::pwd_from_stdin::genome::*;
use crate::pwd_from_stdin::parser;
use crate::pwd_from_stdin::logger;

use crate::pwd_from_stdin::io::*;

use clap::Parser;

use std::fs;
use std::error::Error;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::{self, BufReader, BufRead};
use std::process;

#[macro_use]
extern crate log;


/// Generate a vector of structs `Comparison`. This is the main function enabling batch mode.
/// 
/// TODO: - Convert this into a struct and implement `Display`.
///       - Put this into pileup.rs, or create a separate comparison.rs library.
fn parse_comparisons<'a>(individuals: &[usize], min_depths: &'a [u16], names: &'a [String], allow_self_comparison: bool, genome: &[Chromosome], blocksize: u32) -> Vec<pileup::Comparison> {

    let mut inds = vec![];
    for (i, index) in individuals.iter().enumerate() {
        let name = names.get(i);
        let min_depth = min_depths[(i % (min_depths.len())) as usize]; // wrap around min_depths if its length is lesser than the number of inds.
        inds.push(pileup::Individual::new(name, index, &min_depth));
    }

    let mut comparisons: Vec<pileup::Comparison> = vec![];
    for pair in inds.iter().combinations_with_replacement(2) {
        let self_comparison = pair[0] == pair[1]; 
        if self_comparison && !allow_self_comparison {
            continue
        }
    comparisons.push(pileup::Comparison::new((pair[0].clone(), pair[1].clone()), self_comparison, genome, blocksize));
    }
    comparisons
}

// Main function. 
fn run(cli: parser::Cli) -> Result<(), Box<dyn Error>>{

    // Create output workspace and obtain predefined files.
    let pwd_prefix = cli.get_results_file_prefix()?;

    // ----------------------------- Parse Requested_samples
    let requested_samples: Vec<usize> = parser::parse_user_ranges(&cli.samples, "samples")?;

    // ----------------------------- Sanity checks.
    cli.check_depth()?; // Ensure min_depths are > 2 when allowing self-comparisons
    cli.check_input()?; // Ensure the user has either requested stdin or --pileup

    if cli.min_depth.len() < requested_samples.len() {
        warn!("--min-depth length is less than that of --samples. Values of min-depth will wrap around.")
    }

    // ----------------------------- Initialize genome.
    info!("Indexing reference genome...");
    let genome = match &cli.genome {
        Some(file) => fasta_index_reader(file)?,
        None => default_genome(),
    };

    // ----------------------------- Parse Comparisons
    info!("Parsing Requested comparisons...");
    let mut comparisons = parse_comparisons(&requested_samples, &cli.min_depth, &cli.sample_names, cli.self_comparison, &genome, cli.blocksize);
    let block_prefixes = cli.get_blocks_output_files(&comparisons)?;

    // ----------------------------- Parse target_positions
    info!("Parsing target_positions...");
    let target_positions = match cli.targets {
        None => HashSet::new(),
        Some(filename) => SNPReader::new(&filename)?.hash_target_positions()?
    };

    let target_required: bool = ! target_positions.is_empty();

    // ----------------------------- Parse requested Chromosomes
    let valid_chromosomes : Vec<u8> = match &cli.chr {
        None         => genome.into_iter().map(|chr| chr.name).collect(),
        Some(vector) => parser::parse_user_ranges(vector, "chr")?
    };
    info!("Valid chromosomes: {:?}", valid_chromosomes);

    // ---------------------------- Choose between file handle or standard input
    info!("Opening pileup...");   
    let pileup_reader: Box<dyn BufRead> = match &cli.pileup {
        None => Box::new(BufReader::new(io::stdin())),
        Some(filename) => Box::new(BufReader::new(fs::File::open(filename).unwrap()))
    };

    // ---------------------------- Read Pileup
    info!("Parsing pileup...");   
    for entry in pileup_reader.lines() {
        // ----------------------- Parse line.
        let mut line: pileup::Line = match pileup::Line::new(entry.as_ref().unwrap(), cli.ignore_dels){
            Ok(line) => line,
            Err(e) => {error!("Error: {}", e); process::exit(1);},
        };
        trace!("{:?}", &line);

        // ----------------------- Check if line should be skipped.
        if ! valid_chromosomes.contains(&line.coordinate.chromosome) {
            continue  // Skip line if this is not a valid chromosome. 
        }

        if target_required && !target_positions.contains(&line.coordinate) {  
            continue  // Skip line if we're targeting snps + the current line is not found. 
        }

        // ------------------------ Apply quality filtering on all individuals.
        line.filter_base_quality(&cli.min_qual);

        // ------------------------ Apply target filtration if requested.
        if cli.known_variants {
            let current_coord = match target_positions.get(&line.coordinate) {
                Some(coordinate) => coordinate,
                None => panic!("Cannot filter known variants when REF/ALT allele are unknown! Please use a different file format."),
            };
            line.filter_known_variants(current_coord);
        }

        // ----------------------- Compute PWD (or simply print the line if there's an existing overlap)        
        let mut print_site: bool = false;
        for comparison in &mut comparisons {
            if comparison.satisfiable_depth(&line.individuals) {
                if ! cli.filter_sites {
                    comparison.compare(&line);
                } else {
                    print_site = true;
                }
                if cli.filter_sites &&  print_site {
                    println!("{}", entry.as_ref().unwrap());
                }
            }
        }
    }
    

    // ----------------------------- Print results
    let mut pwd_writer = Writer::new(Some(pwd_prefix["pwd"].clone()))?;
    if ! cli.filter_sites {
        info!("Printing results...");
        let header = format!("{: <20} - Overlap - Sum PWD - Avg. Pwd - Avg. Phred", "Name");
        println!("{}", header);
        pwd_writer.write_iter(&vec![header])?;  // Print to file.
        pwd_writer.write_iter(&comparisons)?;   // 
        for comparison in &comparisons {
            println!("{}", comparison);
        }

        if cli.print_blocks {
            for comparison in &comparisons {
                let pair = comparison.get_pair();
                let mut block_writer = Writer::new(Some(block_prefixes[&pair].clone()))?;
                block_writer.write_iter(vec![&comparison.blocks])?;
                //println!("{}", comparison.blocks);
            }
        }


    }
    Ok(())
}


/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();
    // ----------------------------- Init logger.
    logger::init_logger(&(cli.verbose+(!cli.quiet as u8)));

    // ----------------------------- Serialize command line arguments
    cli.serialize();

    // ----------------------------- Run PWD_from_stdin.
    match run(cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    fn factorial(n: i32 ) -> i32 {
        match n {
            0 => 1,
            1.. => (1..n+1).product(),
            _ => panic!("Not an unsigned number")
        }
    }
    fn permutations_sample_2(n: i32) -> i32 {
        match n {
            0.. => (factorial(n)/(factorial((n-2).abs())*2)),
            _ => panic!()
        }
    }

    #[test]
    fn comparison_length_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let min_depths = vec![2];
            let names = vec![];
            let comparisons = parse_comparisons(&ind_set, &min_depths, &names, true, &default_genome(), 50_000_000);
            let len = ind_set.len() as i32;
            println!("{:?}", comparisons);
            assert_eq!(comparisons.len() as i32, permutations_sample_2(len)+len); 
        }
    }

    #[test]
    fn comparison_length_no_self_allowed() {
        for ind_set in (1..10).powerset().collect::<Vec<_>>().iter() {
            let min_depths = vec![2];
            let names = vec![];
            let comparisons = parse_comparisons(&ind_set, &min_depths, &names, false, &default_genome(), 50_000_000);
            let len = ind_set.len() as i32;
            assert_eq!(comparisons.len() as i32, permutations_sample_2(len)); 
        }
    }
}
