use std::process;
use clap::Parser;

#[macro_use]
extern crate log;

/// Parse command line arguments and run `pwd_from_stdin::run()`
fn main() {
    // ----------------------------- Run CLI Parser 
    let cli = parser::Cli::parse();
    // ----------------------------- Init logger.
    logger::Logger::init(cli.verbose+(!cli.quiet as u8));
    // ----------------------------- Serialize command line arguments
    cli.serialize();
    
    // ----------------------------- unpack Cli and run the appropriate modules.
    match grups::run(cli) {
        Ok(()) => (),
        Err(e) => {
            error!("{}", e);
            process::exit(1);
        }
    };
}
