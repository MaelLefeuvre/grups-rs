use std::path::PathBuf;
use clap::{Parser, Subcommand, Args};


//use serde_yaml;
use serde::{Serialize};

use log::{info};

use std::str::FromStr;
use std::path::{Path};
use core::ops::{Add, Range};
use num::One;

#[derive(Debug)]
pub enum ParserError<'a> {
    ParseIntError(&'a str, String),
    InsufficientDepthError,
    MissingPileupInput
}
impl<'a> std::error::Error for ParserError<'a> {}

impl<'a> std::fmt::Display for ParserError<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::ParseIntError(arg, err) => write!(f, "Invalid slice or integer format for --{}. [{}]", arg, err),
            Self::InsufficientDepthError  => write!(f,"--min_depth must be greater than 1 when performing self-comparison "),
            Self::MissingPileupInput      => write!(f,"Neither --pileup, nor the stdin buffer are being sollicited."),
        }
    }
}

#[derive(Parser, Debug, Serialize)]
pub struct Cli {
    ///Set the verbosity level (-v -vv -vvv -vvvv)
    /// 
    /// Set the verbosity level of this program. With multiple levels
    ///    -v : Info  |  -vv : Debug  | -vvv : Trace
    /// By default, the program will still output Warnings. Use --quiet/-q to disable them
    #[clap(short='v', long, parse(from_occurrences), global=true)]
    pub verbose: u8,
    /// Disable warnings.
    /// 
    /// By default, warnings are emmited and redirected to the console, even without verbose mode on.
    /// Use this argument to disable this. Only errors will be displayed.
    #[clap(short='q', long, global=true)]
    pub quiet: bool,
    #[clap(subcommand)]
    pub commands: Commands,
}

impl<'a> Cli{
    pub fn serialize(&self){
        let serialized = serde_yaml::to_string(&self).unwrap();
        info!("\n---- Command line args ----\n{}\n---", serialized);
    }
}

#[derive(Subcommand, Debug, Serialize)]
pub enum Commands {
    Run {
        #[clap(flatten)]
        common: Common,
        #[clap(flatten)]
        pwd: PwdFromStdin,
        #[clap(flatten)]
        ped: PedigreeSims

    },
    PwdFromStdin {
        #[clap(flatten)]
        common: Common,
        #[clap(flatten)]
        pwd: PwdFromStdin
    },
    PedigreeSims {
        #[clap(flatten)]
        common: Common,
        #[clap(flatten)]
        ped: PedigreeSims
    },
    FST {
        #[clap(flatten)]
        fst: Fst
    }
}

#[derive(Parser, Debug, Serialize)]
pub struct Common {
    /// Minimal required Base Quality (BQ) to perform comparison.
    /// 
    /// Nucleotides whose BQ is lower than the provided treshold are filtered-out. 
    /// Value should be expressed in scale PHRED-33.
    #[clap(short('M'), long, default_value("30"))]
    pub min_qual: u8,
    /// Fasta indexed reference genome.
    /// 
    /// Note that a '.fasta.fai' genome index file must be present at the same directory.
    #[clap(short, long, required(false))]
    pub genome: Option<String>,
    /// Provide with a list of SNP coordinates to target within the pileup.
    /// 
    /// Pileup positions which are not found within the provided --targets file will be exluded from comparison.
    /// Accepted file formats:
    ///   '.snp' : EIGENSTRAT format: space-separated
    ///   '.vcf' : Variant Call Format: tab-separated  
    ///   '.tsv' :   tab-separated file with columns 'CHR POS REF ALT'  
    ///   '.csv' : comma-separated file with columns 'CHR POS REF ALT'  
    ///   '.txt' : space-separated file with columns 'CHR POS REF ALT'   
    #[clap(short, long, required(false))]
    pub targets: Option<String>,
    /// Input pileup file.
    /// 
    /// Note that in the absence of a '--pileup' argument, the program will except a data stream from the standard input. i.e
    /// 'samtools mpileup -B -q30 -Q30 ./samples/* | pwd_from_stdin [...]
    #[clap(short, long, required(false))]
    pub pileup: Option<String>,
    /// Restrict comparison to a given set of chromosomes.
    /// 
    /// Argument may accept slices (inclusive) and/or discrete integers i.e.
    /// '--chr 9-11 13 19-22 ' will be parsed as: [9,10,11,13,19,20,21,22]
    #[clap(short, long, multiple_values(true))]
    pub chr: Option<Vec<String>>,
    /// Provide with a list of sample names for printing.
    /// 
    /// By default, individuals will be referred to using their pileup index. e.g. "Ind0, Ind1, Ind2, etc."
    /// 
    /// The provided list must of course follow the sorted index order which was provided by '--samples'.
    /// Note that the length of --samples and --sample-names do not need to match. When such is the case,
    /// default values are used for the missing names.
    /// 
    /// Example: '--sample 0-2 5 6 --sample-names ART04 ART16 ART20'
    ///   - index: 0      1      2      5     6
    ///   - name : ART04  ART16  ART20  Ind5  Ind6
    #[clap(short='n', long, multiple_values(true))]
    pub sample_names: Vec<String>,
    // Output directory where results will be written.
    #[clap(short, long, default_value("grups-output"))]
    pub output_dir: String,
    //Overwrite existing output files.
    #[clap(short='w', long)]
    pub overwrite: bool,
}

#[derive(Parser, Debug, Serialize)]
pub struct PwdFromStdin {
    /// Enable self-comparison mode on individuals.
    /// 
    /// Note that --min-depth>=2 is required when comparing individuals to themselves, since at least two
    /// alleles are required to effectively compute pairwise differences at a given site.
    #[clap(short('S'), long)]
    pub self_comparison: bool,
    /// Filter out tri-allelic sites when given a list of SNP-coordinates.
    /// 
    /// Note that this arguement requires the use of --targets to provide the program with a list of coordinates.
    /// The given coordinate file must furthermore explicitely provide with the REF/ALT alleles.
    #[clap(short, long)]
    pub known_variants: bool,
    /// Do not perform comparison, but rather print out the pileup lines where a comparison could be made.
    ///
    /// Note that lines are printed as soon as there is a valid comparison between ANY pair of individuals.
    /// Thus, may not prove very effective when using --self-comparison.
    #[clap(short, long)]
    pub filter_sites: bool,
    /// Ignore deletions when performing comparison.
    /// 
    /// By default, deletions ('*' character) will count as a selectible nucleotide for comparison.
    #[clap(short, long)]
    pub ignore_dels: bool,
    /// Print jackknife blocks for each individual
    #[clap(short='J', long)]
    pub print_blocks: bool,
    /// Change the size of our jackknife blocks.
    #[clap(short, long, required(false), default_value("1000000"))]
    pub blocksize: u32,
    
    /// Provide with the minimal sequencing depth required to perform comparison.
    /// 
    /// Note that the length of the provided vector does not need to be the same as the one
    /// provided for individuals. When such is the case, values of --min-depth will "wrap" 
    /// around, and start again from the beginning for the next individuals.
    /// Example:  '--samples 0-4 --min-depth 2 3 5 ' will result in:
    ///   - Ind  : 0 1 2 3 4
    ///   - Depth: 2 3 5 2 3 
    #[clap(short='x', long, multiple_values(true), required(false), default_values(&["1","1"]))]
    pub min_depth: Vec <u16>,
    /// 0-based column index of the individuals that should be compared
    /// 
    /// Argument may accept slices and/or discrete integers i.e.
    /// '--samples 1-4 7 8' will be parsed as: [1, 2, 3, 4, 7, 8]
    /// Note that slices are inclusive.
    #[clap(short, long, multiple_values(true), default_values(&["0", "1"]))]
    pub samples: Vec<String>,
}

#[derive(Args, Debug, Serialize)]
pub struct PedigreeSims {
    /// Proportion of SNPs to keep at true frequency
    #[clap(short='d', long, default_value("0.0"))]
    pub af_downsampling_rate: f64,

    /// Proportion of filtered SNP positions to include in the analysis
    #[clap(short='D', long, default_value("0.0"))]
    pub snp_downsampling_rate: f64,

    /// Contamination rates (or rate ranges) for each desired scenario. 
    /// 
    /// Format: FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ...
    #[clap(short='Q', long, required(false), multiple_values(true), default_values(&["0", "0"]), parse(try_from_str=parse_user_range))]
    pub contamination_rate: Vec<Vec<u8>>,

    /// Sequencing error rates  (or rate ranges) for each desired scenario.
    /// 
    /// Format: FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ...
    /// Default_values: ['FromPileup', 'FromPileup']
    #[clap(short='U', long, required(false), multiple_values(true), default_values(&["0", "0"]), parse(try_from_str=parse_user_range))]
    pub pmd_rate : Vec<Vec<u8>>,

    /// Mean sequencing depths (or depth ranges) for each desired scenario. 
    /// 
    /// Format         : FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ..."
    /// Default values : 'FromPileup'
    #[clap(short='X', long, multiple_values(true), parse(try_from_str=parse_user_range))]
    pub depth_rate: Option<Vec<Vec<u8>>>,

    /// Number of replicates to perform when randomly drawing values from specified ranges for parameters c_rate, mean_cov, q_rate.
    /// 
    /// Format: INT or INT INT INT ... 
    #[clap(short='r', long, required(false), multiple_values(true), default_values(&["1"]))]
    pub param_num_rep: Vec<u32>,
    /// Data output labels.
    /// 
    /// Format: String String String
    #[clap(short, long, required(false), multiple_values(true))]
    pub labels: Vec<String>,

    /// Path to input VCF genomes for founder individuals (1000G)
    #[clap(short='F', long, default_value(r#"./data/founders/1000G-phase3-v5a"#), parse(from_os_str))]
    pub data_dir: std::path::PathBuf,

    /// Path to recombination genetic map data.
    #[clap(short='G', long, required(false), default_value("./data/recombination/GRCh37"), parse(from_os_str))]
    pub recomb_dir: std::path::PathBuf,

    /// Number of pedigree replicates to perform
    #[clap(short='R', long, default_value("1"))]
    pub reps: u32,

    /// Minimum allele frequency within the pedigree superpopulation that is required to include a given SNP within the simulations.
    #[clap(short, long, default_value("0.00"))]
    pub maf: f64,

    /// Superpopulation with which to perform pedigree simulations
    #[clap(short='P', long, default_value("EUR"))]
    pub pedigree_pop: String,

    /// Superpopulation with which to contaminate pedigree simulations.
    /// 
    /// Format: STR, STR
    #[clap(short='C', long, multiple_values(true), max_values(2), default_values(&["AFR","AFR"]))]
    pub contam_pop: Vec<String>,

    ///Number of random individual genomes with which to contaminate pedigree simulations.
    /// 
    /// Format: INT, INT 
    /// Default: Contaminate with population allele frequencies.
    #[clap(short='N', long, multiple_values(true), max_values(2), default_values(&["1","1"]))]
    pub contam_num_ind: Vec<u32>,

    /// Path to input pedigree definition file.
    #[clap(short='T', long, required(false), default_value("./data/pedigrees/default.ped"), parse(from_os_str))]
    pub pedigree: std::path::PathBuf,
    
    // Path to input Panel Definition files.
    #[clap(long, parse(from_os_str))]
    pub panel: Option<std::path::PathBuf>,

    /// Number of parallel CPU processes when performing pedigree-simulations.
    /// 
    /// Parallelization is dispatched according to the number of replicates.
    #[clap(short='@', long, default_value("1"))]
    pub threads: usize,

    ///Number of additional parallel decompression threads when decompressing (BGZF compressed files only).
    /// 
    /// Parallelization is dispatched according to the number of replicates.
    #[clap(short='#', long, default_value("0"))]
    pub decompression_threads: usize

}

#[derive(Args, Debug, Serialize)]
pub struct Fst {

    /// Path to input Panel Definition files.
    #[clap(long, parse(from_os_str))]
    pub panel: std::path::PathBuf,

    /// Population Subset
    /// 
    /// Subset the index by a given number of (super)-population. (e.g. EUR, AFR). Note that at least one pedigree population and one contaminating population
    /// are required to obtain valid index files.
    #[clap(long, multiple_values(true))]
    pub pop_subset: Option<Vec<String>>,

    /// Output directory
    #[clap(long, parse(from_os_str))]
    pub output_dir: std::path::PathBuf,



    /// Path to input VCF genomes for founder individuals (1000G)
    #[clap(short='F', long, default_value(r#"./data/founders/1000G-phase3-v5a"#), parse(from_os_str))]
    pub vcf_dir: std::path::PathBuf,

    /// Number of parallel CPU processes when performing FST-Indexation
    /// 
    /// Parallelization is dispatched according to the number of separate vcf(.gz) files.
    #[clap(short='@', long, default_value("1"))]
    pub threads: usize,

    ///Number of additional parallel decompression threads when decompressing (BGZF compressed files only).
    /// 
    /// Parallelization is dispatched according to the number of replicates.
    #[clap(short='#', long, default_value("0"))]
    pub decompression_threads: usize

}

impl PwdFromStdin {
    /// Sanity check : depth must indeed be > 2 when performing self-comparison.
    /// TODO: - put this in Comparison::new() ?? -> This would allow mix-matching batch mode and self-comparison,
    ///         but could be a bit confusing for users..
    pub fn check_depth<'a>(&self) -> Result<(),ParserError<'a>> {
        if self.self_comparison && self.min_depth.iter().any(|&x| x < 2) {        
            return Err(ParserError::InsufficientDepthError)
        }
        Ok(())
    }
}

/// Compute the average PairWise Difference (PWD) within a pileup file.


/// Command line interface argument parser.
/// 
/// TODO : - `get_results_file_prefix()` and `get_blocks_output_files()` should not be the responsability
///          of this struct. --> migrate to `pwd_from_stdin::io.rs`
///        - add deserialization method. Users could thus fully reproduce a previous run with ease. keep it FAIR. 
impl Common {
    /// Sanity Check: The program should leave if the user did not provide any pileup input, either through
    /// `--pileup` or through stdinput. Without this, our program would wait indefinitely for the stdin buffer.
    pub fn check_input<'a>(&self) -> Result<(), ParserError<'a>> {
        if atty::is(atty::Stream::Stdin) && self.pileup == None {
            return Err(ParserError::MissingPileupInput)
        }
        Ok(())
    }

    /// Get a generic filename for our output files. If the user used `--pileup`, this will become its file stem.
    /// If the user used stdin, this will become a generic name -> "pwd_from_stdin-output" 
    ///
    /// # TODO: This function should be the one responsible of defining the default filename. Stay dry.
    pub fn get_file_prefix(&self, subdir: Option<&str>) -> Option<PathBuf> {
        let default_prefix = String::from("pwd_from_stdin-output");
        let file_prefix = Path::new(self.pileup.as_ref().unwrap_or(&default_prefix)).file_stem()?;

        // Get the path/filename-prefix of all of our outputs.
        let mut parsed_file = PathBuf::new();
        parsed_file.push(&self.output_dir);
        parsed_file.push(subdir.unwrap_or(""));
        parsed_file.push(file_prefix);

        Some(parsed_file)
    }

    /// Check if a given file already exists ; raise an error if such is the case, and the user did not explicitly 
    /// allow file overwriting.
    pub fn can_write_file(&self, pathbuf: &Path) -> std::io::Result<bool> {
        if ! self.overwrite && pathbuf.exists() {   // Check if this file already exists and/or if overwrite is allowed.
            return Err(std::io::Error::new(std::io::ErrorKind::AlreadyExists,
                format!("{:?} exists. use --overwrite to force.",
                pathbuf.to_str().unwrap()))
            )
        }
        Ok(true)
    }
}

/// Convert a user-defined string "range" into a vector of integers.
/// "9-14" thus becomes [9, 10, 11, 12, 13, 14]
/// Note that the range is fully inclusive. 
fn parse_user_range<T>(s: &str) -> Result<Vec<T>, <T as FromStr>::Err> 
where
    T: FromStr + Add<Output = T> + Ord + One,
    Range<T>: Iterator<Item = T>,
{
    match s.split_once('-') {
            Some(t) => Ok((t.0.parse::<T>()?..t.1.parse::<T>()?+One::one()).collect::<Vec<T>>()),
            None    => Ok(vec![s.parse::<T>()?])
    }
}

/// Convert a vector of Strings with user-input ranges to a single vector of integers.
///  - Input will most likely stem from the command line parser, where users are not expected
///    to write every single value they would like to input.
/// 
///     --> ["1-6", "8"] for the user, becomes [1, 2, 3, 4, 5, 6, 8] for our program.
/// 
///  - Return value: Result <
///        - Vector of generic integers, boxed within a Result.
///        - ParseIntError
///    >
/// 
/// # TODO:
///   - Please make this more elegant? We could use a flat_map at some point...
/// 
/// # Example
///```
/////use pwd_from_stdin::parser::parse_user_ranges;
/////let user_input:   Vec<&str> = vec!["5", "1-3", "7"];
/////let parsed_input: Vec<usize>  = parse_user_ranges(&user_input, "chr").unwrap();
/////assert_eq!(parsed_input, vec![1, 2, 3, 5, 7])
///```
/// 
pub fn parse_user_ranges<T>(ranges: Vec<String>, arg_name: &str) -> Result<Vec<T>, ParserError>
where
    T: FromStr + Add<Output = T> + Ord + One,
    <T as FromStr>::Err: ToString,
    Range<T>: Iterator<Item = T>,
{
    let parsed_ranges = ranges.iter().map(|s| parse_user_range(s))
        .collect::<Result<Vec<Vec<T>>, _>>();

    // FIXME - This is inefficient: We're  performing a two-pass iteration while we could 
    //         flatten the structure immediately. use .flat_map() 
    let mut parsed_ranges: Vec<T> = match parsed_ranges {
        Ok(vec) => vec.into_iter().flatten().collect(),
        Err(err) => return Err(ParserError::ParseIntError(arg_name, err.to_string())),
    };

    parsed_ranges.sort();
    parsed_ranges.dedup();
    Ok(parsed_ranges)
}