use clap::{Parser, Subcommand, Args};


#[derive(Parser, Debug)]
struct Cli {
    #[clap(subcommand)]
    commands: Commands,
}

#[derive(Args, Debug)]
struct Common {
    ///Set the verbosity level (-v -vv -vvv -vvvv)
    /// 
    /// Set the verbosity level of this program. With multiple levels
    ///    -v : Info  |  -vv : Debug  | -vvv : Trace
    /// By default, the program will still output Warnings. Use --quiet/-q to disable them
    #[clap(short, long, parse(from_occurrences))]
    pub verbose: u8,
    /// Disable warnings.
    /// 
    /// By default, warnings are emmited and redirected to the console, even without verbose mode on.
    /// Use this argument to disable this. Only errors will be displayed.
    #[clap(short, long)]
    pub quiet: bool,
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
    #[clap(short, long, multiple_values(true))]
    pub sample_names: Vec<String>,
    // Output directory where results will be written.
    #[clap(short, long, default_value("grups-output"))]
    pub output_dir: String,
    //Overwrite existing output files.
    #[clap(short='x', long)]
    pub overwrite: bool,
}

#[derive(Args, Debug)]
struct PwdFromStdin {
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
    #[clap(short='B', long)]
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
    #[clap(short, long, multiple_values(true), required(false), default_values(&["1","1"]))]
    pub min_depth: Vec <u16>,
}

#[derive(Args, Debug)]
struct PedigreeSims {
    #[clap(short, long, required(false))]
    reps: Option<u32>
}

#[derive(Subcommand, Debug)]
enum Commands {
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
    }
}


fn main() {
    println!("Hello, world!");

    let cli = Cli::parse();


    match &cli.commands {
        Commands::Run {common:_, pwd:_, ped} => {
            println!("Command Run : {:#?}", cli.commands);
            println!("Hi!  : {:#?}", ped.reps);

        },
        Commands::PwdFromStdin {common:_, pwd} => {
            println!("Command pwd_from_stdin: {:#?}", pwd.min_depth);
        },
        Commands::PedigreeSims {..} => {
            println!("Command pedigree-sims: {:#?}", cli.commands);
        }
    };
}
