extern crate parser;
extern crate logger;

use parser::{Cli, Commands::*};
use genome::Genome;

#[macro_use]
extern crate log;

use std::error::Error;

/// @TODO : Stay dry...
pub fn run(cli: Cli) -> Result<(), Box<dyn Error>> {
    match cli.commands {
        PedigreeSims {common, pwd, ped} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None => Genome::default(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (mut comparisons, _target_positions) = pwd_from_stdin::run(&common, &pwd, &requested_samples, &genome)?;
            let _ = pedigree_sims::run(common, ped, &mut comparisons)?;
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
        // PedigreeSims {..} => {
        //    println!("Command pedigree-sims: {:#?}", cli.commands);
        //},

        FST {fst: fst_cli} => {
            vcf_fst::run(&fst_cli)?
        }
    };
    Ok(())
}