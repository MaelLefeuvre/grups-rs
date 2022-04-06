use clap::Parser;


#[derive(Parser, Debug)]
pub struct Cli {
    #[clap(short, long)]
    pub verbose: bool,
    #[clap(short, long)]
    pub quiet: bool,
    #[clap(short('S'), long)]
    pub self_comparison: bool,
    #[clap(short, long)]
    pub known_variants: bool,
    //#[clap(short, long)]
    //pub filter_sites: bool,
    #[clap(short, long)]
    pub ignore_dels: bool,
    #[clap(short, long)]
    pub print_blocks: bool,
    #[clap(short, long, required(false), default_value("1000000"))]
    pub blocksize: u32,
    //#[clap(short, long, number_of_values(2), required(false), default_values(&["1","1"]))]
    #[clap(short, long, multiple_values(true), required(false), default_values(&["1","1"]))]
    pub min_depth: Vec <u16>,
    #[clap(short('M'), long, default_value("30"))]
    pub min_qual: u8,
    #[clap(short, long, multiple_values(true), default_values(&["0", "1"]))]
    pub samples: Vec <usize>,
    //#[clap(short, long, required(false), default_value("1000000"))]
    //pub lines: u32,
    #[clap(short, long, required(false))]
    pub targets: Option<String>,
    #[clap(short, long, required(false))]
    pub pileup: Option<String>,
    #[clap(short, long, multiple_values(true))]
    pub chr: Option<Vec<u8>>,
    #[clap(short, long, multiple_values(true))]
    pub sample_names: Vec<String>,
}
