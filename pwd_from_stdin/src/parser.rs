use clap::Parser;

extern crate serde_yaml;
use serde::{Serialize};

use log::{info};



#[derive(Parser, Debug, Serialize)]
pub struct Cli {
    #[clap(short, long, parse(from_occurrences))]
    pub verbose: u8,
    #[clap(short, long)]
    pub quiet: bool,
    #[clap(short('S'), long)]
    pub self_comparison: bool,
    #[clap(short, long)]
    pub known_variants: bool,
    #[clap(short, long)]
    pub filter_sites: bool,
    #[clap(short, long)]
    pub ignore_dels: bool,
    #[clap(short, long)]
    pub print_blocks: bool,
    #[clap(short, long, required(false), default_value("1000000"))]
    pub blocksize: u32,
    #[clap(short, long, multiple_values(true), required(false), default_values(&["1","1"]))]
    pub min_depth: Vec <u16>,
    #[clap(short('M'), long, default_value("30"))]
    pub min_qual: u8,
    #[clap(short, long, multiple_values(true), default_values(&["0", "1"]))]
    pub samples: Vec <usize>,
    //#[clap(short, long, required(false), default_value("1000000"))]
    //pub lines: u32,
    #[clap(short, long, required(false))]
    pub genome: Option<String>,
    #[clap(short, long, required(false))]
    pub targets: Option<String>,
    #[clap(short, long, required(false))]
    pub pileup: Option<String>,
    #[clap(short, long, multiple_values(true))]
    pub chr: Option<Vec<u8>>,
    #[clap(short, long, multiple_values(true))]
    pub sample_names: Vec<String>,
}

impl Cli {
    pub fn serialize(&self){
        let serialized = self::serde_yaml::to_string(&self).unwrap();
        info!("\n---- Command line args ----\n{}\n---", serialized);
    }

}
