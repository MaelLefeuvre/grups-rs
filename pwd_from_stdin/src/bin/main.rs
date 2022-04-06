extern crate pwd_from_stdin;

use crate::pwd_from_stdin::pileup;
//use crate::pwd_from_stdin::jackknife::*;
use crate::pwd_from_stdin::genome::*;
use crate::pwd_from_stdin::parser::Cli;


use clap::Parser;


use std::{env, fs};
use std::error::Error;
use itertools::Itertools;
use std::collections::HashSet;
use std::io::{self, BufReader, BufRead};


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


fn default_genome() -> Vec<Chromosome> {
    vec![
        Chromosome{index:  0, name:  1, length: 249250621},
        Chromosome{index:  1, name:  2, length: 243199373},
        Chromosome{index:  2, name:  3, length: 198022430},
        Chromosome{index:  3, name:  4, length: 191154276},
        Chromosome{index:  4, name:  5, length: 180915260},
        Chromosome{index:  5, name:  6, length: 171115067},
        Chromosome{index:  6, name:  7, length: 159138663},
        Chromosome{index:  7, name:  8, length: 146364022},
        Chromosome{index:  8, name:  9, length: 141213431},
        Chromosome{index:  9, name: 10, length: 135534747},
        Chromosome{index: 10, name: 11, length: 135006516},
        Chromosome{index: 11, name: 12, length: 133851895},
        Chromosome{index: 12, name: 13, length: 115169878},
        Chromosome{index: 13, name: 14, length: 107349540},
        Chromosome{index: 14, name: 15, length: 102531392},
        Chromosome{index: 15, name: 16, length:  90354753},
        Chromosome{index: 16, name: 17, length:  81195210},
        Chromosome{index: 17, name: 18, length:  78077248},
        Chromosome{index: 18, name: 19, length:  59128983},
        Chromosome{index: 19, name: 20, length:  63025520},
        Chromosome{index: 20, name: 21, length:  48129895},
        Chromosome{index: 21, name: 22, length:  51304566}
    ]
}



fn main() {
    let cli = Cli::parse();
    // ----------------------------- Initialize defaults
    let genome = default_genome();

    // ----------------------------- Command line arguments   --> This is horrible. put a method in
    //                                                            parser.rs
    //println!("Filter sites: {}", cli.filter_sites);
    println!("requested_individuals: {:?}", cli.samples);
    println!("requested_min_depth: {:?}", cli.min_depth);
    println!("Allow self comparison: {}", cli.self_comparison);
    println!("Phred treshold: {}", cli.min_qual);
    println!("ignore_dels: {:?}", cli.ignore_dels);
    println!("Filter known_variants: {}", cli.known_variants);
    println!("Print Jackknife blocks: {}", cli.print_blocks);
    println!("Jackknife blocksize: {}", cli.blocksize);
    println!("Input file: {:?}", cli.pileup.as_ref());
    println!("Targets file: {:?}", cli.targets.as_ref());

    let sep = " ";

    // ----------------------------- Sanity checks!
    if cli.self_comparison && cli.min_depth.iter().any(|&x| x < 2) {         // depth must be > 2 when performing self-comparison
        panic!("Min_depth must be greater than 1 when performing self-comparison");
    }

    // ----------------------------- Parse jackknife blocks
    //let jackknife_blocks = compute_jackknife_blocksize(&genome, blocksize);
    //let jackknife_blocks = JackknifeBlocks::new(&genome, blocksize);

    // ----------------------------- Parse Comparisons
    //let requested_names = vec!["MT23", "MT26", "MT7"];
    let mut comparisons = parse_comparisons(&cli.samples, cli.min_depth, cli.sample_names, cli.self_comparison, &genome, cli.blocksize);


    // ----------------------------- Parse target_positions
    println!("// ------------------- Parse target positions -------------------- //");   
    let target_positions = match cli.targets {
        None => HashSet::new(),
        Some(filename) => match hash_target_positions(&filename, &sep) {
            Ok(vector) => vector,
            Err(error) => panic!("Problem parsing the file. {:?}", error),
        },
    };
    let target_required: bool = ! target_positions.is_empty();
    println!("Target_required: {}", target_required);



    // --------------------------- Parse chromosomes 
    println!("// ------------------- Parse Chromosomes      -------------------- //");   
    let valid_chromosomes : Vec<u8> = match cli.chr {
        Some(vector) => vector,
        None => genome.into_iter().map(|chr| chr.name).collect()
    };
    println!("Valid chromosomes: {:?}", valid_chromosomes);
    
    // ---------------------------- Choose between file handle or standard input
    println!("// ------------------- Opening pileup...      -------------------- //");   
    let pileup_reader: Box<dyn BufRead> = match cli.pileup {
        None => Box::new(BufReader::new(io::stdin())),
        Some(filename) => Box::new(BufReader::new(fs::File::open(filename).unwrap()))
    };

    // ---------------------------- Read Pileup
    println!("// ------------------- Parsing pileup...      -------------------- //");   
    for entry in pileup_reader.lines() {
        // ----------------------- Parse line.
        let mut line: pileup::Line = pileup::Line::new(&entry.as_ref().unwrap(), '\t', cli.ignore_dels);

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

        for comparison in &mut comparisons {
            if comparison.satisfiable_depth(&line.individuals) {
                comparison.compare(&line);
                //if filter_sites {
                //    println!("{}", entry.unwrap());
                //} 
            }
        }
    }

    println!("// ------------------- Printing results       -------------------- //");   
    println!("{: <20} - Overlap - Sum PWD - Avg. Pwd - Avg. Phred", "Name");
    for comparison in &comparisons {
        comparison.print();
        if cli.print_blocks {
            comparison.blocks.print();
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
