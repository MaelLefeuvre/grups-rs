//extern crate pwd_from_stdin;

use pwd_from_stdin::parser;
use pwd_from_stdin::logger;

use clap::Parser;

use std::process;

#[macro_use]
extern crate log;

/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();
    // ----------------------------- Init logger.
    logger::init_logger(&(cli.verbose+(!cli.quiet as u8)));

    // ----------------------------- Serialize command line arguments
    cli.serialize();

    // ----------------------------- Run PWD_from_stdin.
    match pwd_from_stdin::run(cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}
