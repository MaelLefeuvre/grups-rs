extern crate pwd_from_stdin;

use crate::pwd_from_stdin::pileup;
//use crate::pwd_from_stdin::jackknife::*;
use crate::pwd_from_stdin::genome::*;
use crate::pwd_from_stdin::parser::Cli;
use crate::pwd_from_stdin::logger;

use clap::Parser;

use std::fs;
use std::error::Error;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::{self, BufReader, BufRead};
use std::process;
use atty;

#[macro_use]
extern crate log;


/// Convert a space-separated path of SNP coordinates to a vector of object SNPCoord.
/// TODO: At this state, does not support multiple spaces. Should implement ability to
///       Remove empty fields.
///
/// Input: - path (string): path and filename to a snp coordinates file.
///          Columns are: CHR (Required)
///                       POS (Required)
///                       REF (Optional)
///                       ALT (Optional)
///
/// Return: Vector of structs 'SNPCoord'
fn hash_target_positions(path: &String, sep: &str) -> Result<HashSet<SNPCoord>, Box<dyn Error>> {
    let mut target_positions : HashSet<SNPCoord> = HashSet::new(); // Output
    let file = BufReader::new(fs::File::open(path).unwrap());
    for line in file.lines() {
        let line = line.unwrap();
        let split_line: Vec<&str>    = line.split(sep).collect();
        let chromosome: u8           = split_line[0].parse().unwrap();
        let position  : u32          = split_line[1].parse().unwrap();
        let reference : Option<char> = split_line[2].parse().ok();
        let alternate : Option<char> = split_line[3].parse().ok();
        
        let coordinate: SNPCoord     = SNPCoord {chromosome, position, reference, alternate};
        target_positions.insert(coordinate);

    }
    Ok(target_positions)
}

fn parse_comparisons<'a>(individuals: &Vec<usize>, min_depths: Vec<u16>, names: Vec<String>, allow_self_comparison: bool, genome: &Vec<Chromosome>, blocksize: u32) -> Vec<pileup::Comparison> {

    let mut inds = vec![];
    for (i, index) in individuals.iter().enumerate() {
        let name = names.get(i);
        let min_depth = min_depths[(i % (min_depths.len())) as usize]; // wrap around min_depths if its length is lesser than the number of inds.
        inds.push(pileup::Individual::new(name, index, &min_depth));
    }

    let mut comparisons: Vec<pileup::Comparison> = vec![];
    for pair in inds.iter().combinations_with_replacement(2) {
        let self_comparison = &pair[0] == &pair[1]; 
        if self_comparison && !allow_self_comparison {
            continue
        }
    comparisons.push(pileup::Comparison::new((pair[0].clone(), pair[1].clone()), self_comparison, genome, blocksize));
    }
    comparisons
}

fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = Cli::parse();
    // ----------------------------- Init logger.
    logger::init_logger(&(cli.verbose+(!cli.quiet as u8)));

    // ----------------------------- Serialize command line arguments
    cli.serialize();

    // ----------------------------- Sanity checks!
    if cli.self_comparison && cli.min_depth.iter().any(|&x| x < 2) {         // depth must be > 2 when performing self-comparison
        error!("Min_depth must be greater than 1 when performing self-comparison");
        process::exit(1);
    }

    if atty::is(atty::Stream::Stdin) && cli.pileup == None {
        error!("Neither --pileup, nor the stdinput buffer are being sollicited. Exiting.");
        process::exit(1);
    }

    if cli.min_depth.len() < cli.samples.len() {
        warn!("--min-depth length is less than that of --samples. Values of min-depth will wrap around.")
    }
    // ----------------------------- Initialize genome.
    info!("Indexing reference genome...");
    let genome = match cli.genome.as_ref(){
        Some(file) => fasta_index_reader(&(file.to_string()+&".fai".to_string())).unwrap(),
        None => default_genome(),
    };

    // ----------------------------- Parse Comparisons
    info!("Parsing Requested comparisons...");
    let mut comparisons = parse_comparisons(&cli.samples, cli.min_depth, cli.sample_names, cli.self_comparison, &genome, cli.blocksize);

    // ----------------------------- Parse target_positions
    info!("Parsing target_positions..."); 
    let sep = " ";
    let target_positions = match cli.targets {
        None => HashSet::new(),
        Some(filename) => match hash_target_positions(&filename, &sep) {
            Ok(vector) => vector,
            Err(error) => panic!("Problem parsing the file. {:?}", error),
        },
    };
    let target_required: bool = ! target_positions.is_empty();

    // --------------------------- Parse chromosomes 
    let valid_chromosomes : Vec<u8> = match cli.chr {
        Some(vector) => vector,
        None => genome.into_iter().map(|chr| chr.name).collect()
    };
    info!("Valid chromosomes: {:?}", valid_chromosomes);
    
    // ---------------------------- Choose between file handle or standard input
    info!("Opening pileup...");   
    let pileup_reader: Box<dyn BufRead> = match cli.pileup {
        None => Box::new(BufReader::new(io::stdin())),
        Some(filename) => Box::new(BufReader::new(fs::File::open(filename).unwrap()))
    };

    // ---------------------------- Read Pileup
    info!("Parsing pileup...");   
    for entry in pileup_reader.lines() {
        // ----------------------- Parse line.
        let mut line: pileup::Line = match pileup::Line::new(&entry.as_ref().unwrap(), '\t', cli.ignore_dels){
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
            line.filter_known_variants(&current_coord);
        }

        //let filter_sites = true;
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
                //comparison.compare(&line);
                //if filter_sites {
                //    println!("{}", entry.unwrap());
                //} 
            }
        }
    }

    if ! cli.filter_sites {
        info!("Printing results...");
        println!("{: <20} - Overlap - Sum PWD - Avg. Pwd - Avg. Phred", "Name");
        for comparison in &comparisons {
            println!("{}", comparison);
            if cli.print_blocks {
                comparison.blocks.print();
            }
        }
    }
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
            let comparisons = parse_comparisons(&ind_set, min_depths, names, true, &default_genome(), 50_000_000);
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
            let comparisons = parse_comparisons(&ind_set, min_depths, names, false, &default_genome(), 50_000_000);
            let len = ind_set.len() as i32;
            assert_eq!(comparisons.len() as i32, permutations_sample_2(len)); 
        }
    }
}
