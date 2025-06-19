use std::{io, fs::File, collections::HashMap, io::{BufReader, BufRead}};

use ahash::AHashSet;

use log::{warn, info, debug};
use located_error::prelude::*;

use genome::{Genome, coordinate::Coordinate};
use grups_io::{parse::{self, FileKey}, read::SNPReader};

pub mod pileup;

pub mod comparisons;
use comparisons::Comparisons;

pub mod error;
pub use error::PwdFromStdinError;

/// Run pwd-from-stdin, using the user-provided parameters and input file.
/// Returns a `Comparisons` struct, containing the results.
/// 
/// # Errors 
/// - If any of the Input/Output file is invalid (`FileNotFound`, `PermissionDenied`, etc.)
/// - If at any point, the program fails to parse a line from the user provided pileup.
pub fn run<'a>(
    com_cli           : &'a parser::Common,
    pwd_cli           : &'a parser::PwdFromStdin,
    requested_samples : &'a [usize],
    genome            : &'a Genome,
) -> Result<(Comparisons, HashMap<String, String>)> {
    info!("Running 'pwd-from-stdin' module...");
    // ----------------------------- Sanity checks.
    pwd_cli.check_depth()?; // Ensure min_depths are > 2 when allowing self-comparisons
    com_cli.check_input()?; // Ensure the user has either requested stdin or --pileup

    if pwd_cli.min_depth.len() < requested_samples.len() {
        warn!("Number of provided --min-depth values is less than that of --samples. Values will be recycled.");
    }

    // ----------------------------- Parse Comparisons
    info!("Parsing Requested comparisons...");
    let mut comparisons = Comparisons::parse(requested_samples, &pwd_cli.min_depth, &com_cli.sample_names, pwd_cli.self_comparison, genome, pwd_cli.blocksize);

    // ----------------------------- Prepare output files
    // ---- Add pwd files.
    let mut output_files = parse::get_output_files(
        &mut com_cli.get_file_prefix(None)?,  // extract the user requested file prefix
        com_cli.overwrite,                    // Should we allow file overwriting ?
        FileKey::Ext,                         // What key are we using to hash these files ?
        &["".to_string()],                    // Vector of filename suffixes.
        &["pwd"]                              // Vector of file extensions.
    )?;

    // ---- Add blocks files.
    output_files.extend(
        parse::get_output_files(
            &mut com_cli.get_file_prefix(Some("blocks/"))?,
            com_cli.overwrite,
            FileKey::Suffix,
            &comparisons.get_pairs(),
            &["blk"]
        )?);

    debug!("Output files: {:#?}", output_files);


    // ----------------------------- Parse target_positions
    info!("Parsing target_positions...");
    let target_positions = match &com_cli.targets {
        None => AHashSet::new(),
        Some(filename) => SNPReader::new(filename)?.hash_target_positions(pwd_cli.exclude_transitions)?
    };

    let target_required: bool = ! target_positions.is_empty();

    // Early exit If the user requested transition filtration without specifying an SNP panel.
    if ! target_required && pwd_cli.exclude_transitions {
        return Err(PwdFromStdinError::MissingTargetPositions).loc("While initializing main event loop")
    }

    // ----------------------------- Parse requested Chromosomes
    let valid_chromosomes : Vec<u8> = match com_cli.chr.clone() {
        None         => genome.keys().copied().map(|x| x.into()).collect(),
        Some(vector) => {
            let vec: Vec<String> = vector.into_iter().map(|c| match c.as_ref() {
                "X" | "chrX" => "88".to_string(),
                other => other.to_string()
            }).collect::<Vec<String>>();
            parser::parse_user_ranges(&vec, "chr")?
        }
    };
    info!("Valid chromosomes: {:?}", valid_chromosomes);

    // ---------------------------- Choose between file handle or standard input
    info!("Opening pileup...");   
    let pileup_reader: Box<dyn BufRead> = match &com_cli.pileup {
        None => Box::new(BufReader::new(io::stdin())),
        Some(filename) => Box::new(BufReader::new(File::open(filename)?))
    };

    // ---------------------------- Read Pileup
    info!("Parsing pileup...");   
    let loc_msg = {|c: &Coordinate| format!("While parsing coordinate coordinate: {c}")};
    for (i, entry) in pileup_reader.lines().enumerate() {
        // ----------------------- Parse line.
        let entry = entry?;
        let mut line = pileup::Line::new(&entry, !pwd_cli.consider_dels)
            .with_loc(|| format!("While attempting to parse pileup line n°{}", i+1))?;


        // ----------------------- Check if line should be skipped.
        if ! valid_chromosomes.contains(&line.coordinate.chromosome) {
            continue  // Skip line if this is not a valid chromosome. 
        }

        if target_required && !target_positions.contains(&line.coordinate) {  
            continue  // Skip line if we're targeting snps + the current line is not found. 
        }

        // ------------------------ Apply quality filtering on all individuals.
        line.filter_base_quality(com_cli.min_qual);

        // ------------------------ Apply target filtration if requested.
        if pwd_cli.known_variants {
            let known_coord = target_positions.get(&line.coordinate)
                .ok_or(PwdFromStdinError::MissingKnownVariant)
                .with_loc(|| loc_msg(&line.coordinate))?;
            line.filter_known_variants(known_coord)?;
        }

        // ----------------------- Compute PWD (or simply print the line if there's an existing overlap)        
        for comparison in comparisons.iter_mut() {
            
            let satisfactory_depth = comparison.satisfiable_depth(&line.individuals)
                .with_context(|| loc_msg(&line.coordinate))
                .with_loc(|| format!("While comparing pair {}", comparison.get_pair()))?;

            if satisfactory_depth {
                if pwd_cli.filter_sites {
                    println!("{}", &entry);
                } else {
                    comparison.compare(&line)?;
                }
            }
        }
    }

    // Run two-pass variance estimation algorithm.
    comparisons.update_variance_unbiased();

    Ok((comparisons, output_files))
}
