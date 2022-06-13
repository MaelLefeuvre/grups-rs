use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();

    // ----------------------------- Init logger.
    let verbosity = if cli.quiet {0} else {cli.verbose+1};
    logger::Logger::init(verbosity);
    
    // ----------------------------- Serialize command line arguments
    match cli.serialize() {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
    
    // ----------------------------- unpack Cli and run the appropriate modules.
    match grups::run(cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}
