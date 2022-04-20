use parser::{Cli, Commands::*};
use pwd_from_stdin::genome::*;

use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

use std::error::Error;


/// @TODO : Stay dry...
fn run(cli: &Cli) -> Result<(), Box<dyn Error>> {
    match &cli.commands {
        Run {common, pwd, ped} => {
            println!("Command Run : {:#?}", cli.commands);
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => fasta_index_reader(file)?,
                None => default_genome(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (comparisons, target_positions) = pwd_from_stdin::run(common, pwd, &requested_samples, &genome)?;
            pedigree_sims::run(common, ped, &requested_samples, &genome, &Some(comparisons), &Some(target_positions))?;
        },
        PwdFromStdin {common, pwd} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => fasta_index_reader(file)?,
                None => default_genome(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (_,_) = pwd_from_stdin::run(common, pwd, &requested_samples, &genome)?;
        },
        PedigreeSims {..} => {
            println!("Command pedigree-sims: {:#?}", cli.commands);
        }
    };
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
    
    // ----------------------------- unpack Cli and run the appropriate modules.
    match run(&cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}
