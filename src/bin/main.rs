use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

use anyhow::Result;

/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() -> Result<()> {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();

    // ----------------------------- Init logger.
    let verbosity = if cli.quiet {0} else {cli.verbose + 1};
    logger::Logger::init(verbosity);
    
    // ----------------------------- Serialize command line arguments
    if let Err(e) = cli.serialize() {
        error!("{:?}", e);
        process::exit(1);
    };
    
    // ----------------------------- unpack Cli and run the appropriate modules.
    if let Err(e) = grups::run(cli) {
        error!("{:?}", e);
        process::exit(1);
    };

    Ok(())
}
