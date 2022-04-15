use logger;
use parser::{Commands::*};

use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();
    // ----------------------------- Serialize command line arguments
    cli.serialize();
    match &cli.commands {
        Run {common:_, pwd:_, ped:_} => {
            println!("Command Run : {:#?}", cli.commands);
        },
        PwdFromStdin {common, pwd} => {
            // ----------------------------- Init logger.
            logger::init_logger(&(common.verbose+(!common.quiet as u8)));

            // ----------------------------- Run PWD_from_stdin.
            match pwd_from_stdin::run(common, pwd) {
                Ok(()) => (),
                Err(e) => {
                    error!("{}", e);
                    process::exit(1);
                }
            };
        },
        PedigreeSims {..} => {
            println!("Command pedigree-sims: {:#?}", cli.commands);
        }
    };
}
