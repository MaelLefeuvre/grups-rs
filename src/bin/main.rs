use parser::{Cli, Commands::*};
use genome::Genome;

use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

use std::error::Error;

/// @TODO : Stay dry...
fn run(cli: Cli) -> Result<(), Box<dyn Error>> {
    match cli.commands {
        Run {common, pwd, ped} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None => Genome::default(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (comparisons, _target_positions) = pwd_from_stdin::run(&common, &pwd, &requested_samples, &genome)?;
            //let comparisons = Arc::new(comparisons);
            pedigree_sims::run(common, ped, genome, &comparisons)?;
        },
        PwdFromStdin {common, pwd} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None => Genome::default(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (_,_) = pwd_from_stdin::run(&common, &pwd, &requested_samples, &genome)?;
        },
        PedigreeSims {..} => {
            println!("Command pedigree-sims: {:#?}", cli.commands);
        },

        FST {fst: fst_cli} => {
            vcf_fst::run(&fst_cli)?
        }
    };
    Ok(())
}


/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();
    // ----------------------------- Init logger.
    logger::Logger::init(cli.verbose+(!cli.quiet as u8));
    // ----------------------------- Serialize command line arguments
    cli.serialize();
    
    // ----------------------------- unpack Cli and run the appropriate modules.
    match run(cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}
