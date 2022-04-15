use logger;
//use parser;
//
//use parser::{PwdFromStdin, Commands};
//
//use std::process;
//use clap::Parser;
//
//#[macro_use]
//extern crate log;
//
///// Parse command line arguments and run `pwd_from_stdin::run()`
//fn main() {
//    // ----------------------------- Run CLI Parser 
//    let cli = parser::PwdFromStdin::parse();
//    // ----------------------------- Init logger.
//    logger::init_logger(&(cli.verbose+(!cli.quiet as u8)));
//
//    // ----------------------------- Serialize command line arguments
//    cli.serialize();
//
//    // ----------------------------- Run PWD_from_stdin.
//    match pwd_from_stdin::run(cli) {
//        Ok(()) => (),
//        Err(e) => {
//            error!("{}", e);
//            process::exit(1);
//        }
//    };
//}
//