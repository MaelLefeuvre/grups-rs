use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

use anyhow::Result;

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;


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
    if let Err(e) = grups_rs::run(cli) {
        error!("{:?}", e);
        process::exit(1);
    };

    Ok(())
}
