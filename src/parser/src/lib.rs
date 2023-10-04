use std::{
    error::Error,
    fs::File,
    path::{Path, PathBuf},
    str::FromStr,
    ops::{Add, Range},
    fmt::{self, Display, Formatter}, ffi::OsStr
};

use located_error::*;

use clap::{Parser, Subcommand, Args, ArgEnum, ArgAction};
use serde::{Serialize, Deserialize};
use log::debug;
use num::One;
use anyhow::{anyhow, Result};

mod error;
pub use error::ParserError;

#[derive(Parser, Debug, Serialize, Deserialize)]
#[clap(name="grups-rs", author, version, about, long_about = None)]
#[clap(propagate_version = true)]
/// GRUPS-rs: Get Relatedness Using Pedigree Simulations
pub struct Cli {
    ///Set the verbosity level (-v -vv -vvv)
    /// 
    /// Set the verbosity level of this program. Multiple levels allowed {n} 
    ///
    /// -v: Info  |  -vv: Debug  | -vvv: Trace {n}
    /// 
    /// Note that the program will still output warnings by default, even when this flag is off.
    /// Use The --quiet/-q to disable them
    #[clap(short='v', long, parse(from_occurrences), global=true)]
    pub verbose: u8,

    /// Disable warnings.
    /// 
    /// By default, warnings are emmited and redirected to the console, even when verbose mode is off.
    /// Use this argument to disable this. Only errors will be displayed.
    #[clap(short='q', long, global=true)]
    pub quiet: bool,

    #[clap(subcommand)]
    pub commands: Commands,
}

impl Cli{
    /// Serialize command line arguments within a `.yaml` file.
    /// 
    /// # Behavior
    /// - File naming follows the convention '{current time}-{module name}.yaml'. current time follows the format
    ///   `YYYY`-`MM`-`DD`T`hhmmss`
    /// - File is written at the root of the user-provided `--output-dir` folder.
    /// 
    /// # Errors
    /// Sends an unrecoverable error if `serde_yaml` fails to parse `Self` to a string.
    /// 
    /// # Panics
    /// - Throws a tantrum whenever the provided `--output-dir` is invalid.
    pub fn serialize(&self) -> Result<(), Box<dyn Error>> {

        // Parse arguments to yaml and print to console.
        let serialized = serde_yaml::to_string(&self)
            .map_err(|err| format!("Failed to serialize command line arguments. got [{err}]"))?;
        
        debug!("\n---- Command line args ----\n{}\n---", serialized);

        // Fetch the appropriate output-directory and parse the name of the output file.
        let current_time = chrono::offset::Local::now().format("%Y-%m-%dT%H%M%S").to_string();

        let output_file = match &self.commands { // @TODO: lots of repetition here. Stay dry: macro_rules!
            Commands::PedigreeSims {common, pwd: _, ped: _} => {
                let dir_string = common.output_dir.to_str().expect("Invalid characters in directory");
                format!("{dir_string}/{current_time}-pedigree-sims.yaml")
            },
            Commands::PwdFromStdin {common, pwd: _} => {
                let dir_string = common.output_dir.to_str().expect("Invalid characters in directory");
                format!("{dir_string}/{current_time}-pwd-from-stdin.yaml")
            },
            Commands::FST {fst} => {
                let dir_string = fst.output_dir.to_str().expect("Invalid characters in directory");
                format!("{dir_string}/{current_time}-fst-index.yaml")
            },

            Commands::FromYaml {yaml: _} => return Ok(()),
            Commands::Cite => return Ok(())
        };

        // Write arguments
        match std::fs::write(&output_file, serialized) {
            Err(e) => Err(format!("Unable to serialize arguments into {output_file}: [{e}]").into()),
            Ok(()) => Ok(()),
        }
    }

    /// Deserialize a `.yaml` file into Command line arguments.
    /// 
    /// # Errors
    /// 
    /// - Returns `FileNotFound` or `PermissionDenied` if the provided `.yaml` is invalid,
    ///   or does not carry read permissions
    /// - Sends an unrecoverable error if: `serde_yaml` fails to parse the provided file to `Self`.
    pub fn deserialize(yaml: PathBuf) -> Result<Self, Box<dyn Error>> {
        Ok(serde_yaml::from_reader(File::open(yaml)?)?)
    }
}

#[derive(Subcommand, Debug, Serialize, Deserialize)]
pub enum Commands {
    PedigreeSims {
        #[clap(flatten)]
        common: Common,
        #[clap(flatten)]
        pwd: PwdFromStdin,
        #[clap(flatten)]
        ped: Box<PedigreeSims> // Box<T> to mitigate the large size difference between variants.

    },
    PwdFromStdin {
        #[clap(flatten)]
        common: Common,
        #[clap(flatten)]
        pwd: PwdFromStdin
    },

    FST {
        #[clap(flatten)]
        fst: VCFFst
    },

    /// Run grups-rs using a previously generated .yaml configuration file.
    /// 
    /// This allows users to easily re-apply a grups-rs command using the exact same parameters
    /// and arguments. 
    FromYaml {
        yaml: PathBuf,
    },

    Cite 
}

/// Print all citations tied to this project
pub struct Cite;

#[derive(Parser, Debug, Default, Serialize, Deserialize)]
pub struct Common {
    /// Minimal required Base Quality (BQ) to perform comparison.
    /// 
    /// Nucleotides whose base quality is lower than the provided treshold will be filtered-out. 
    /// Value should be expressed in Phred-33 scale.
    #[clap(short('M'), long, default_value("30"))]
    pub min_qual: u8,

    /// Fasta indexed reference genome.
    /// 
    /// By default, grups-rs will use GRCh37 as a reference genome. Use this argument if you wish to specify an alternative reference.
    /// Note that a '.fasta.fai' genome index file must be present at the same directory.
    #[clap(short, long, required(false))]
    pub genome: Option<String>,

    /// Provide with a list of SNP coordinates to target within the pileup.
    /// 
    /// Pileup positions which are not found within the provided --targets file will be exluded from comparison.  
    /// 
    /// Accepted file formats:{n}
    ///   '.snp' : EIGENSTRAT format: space-separated{n}
    ///   '.vcf' : Variant Call Format: tab-separated{n}
    ///   '.tsv' :   tab-separated file with columns 'CHR POS REF ALT'{n}
    ///   '.csv' : comma-separated file with columns 'CHR POS REF ALT'{n}
    ///   '.txt' : space-separated file with columns 'CHR POS REF ALT'{n}
    /// 
    #[clap(short, long, required(false))]
    pub targets: Option<String>,

    /// Input Pileup file.
    /// 
    /// Note that in the absence of a '--pileup' argument, the program may accept a data stream from the standard input. i.e:{n} 
    ///
    ///     samtools mpileup -B -q30 -Q30 [samples.bam] | grups-rs pwd_from_stdin [...]{n}
    /// 
    #[clap(long, required(false))]
    pub pileup: Option<String>,

    /// Restrict comparison to a given set of chromosomes.
    /// 
    /// Argument may accept slices (inclusive) such as '--chr 9-11' and/or discrete integers such as '--chr 1 4 13'.{n}
    /// Example:{n}
    ///   specifying          : '--chr 9-11 13 19-22 '{n}
    ///   ...will be parsed as: [9, 10, 11, 13, 19, 20, 21, 22]
    /// 
    #[clap(short, long, multiple_values(true))]
    pub chr: Option<Vec<String>>,

    /// Provide with a list of sample names for printing.{n}
    /// 
    /// By default, individuals will be referred to using their pileup index. e.g. "Ind0, Ind1, Ind2, etc."
    /// 
    /// The provided list must of course follow the sorted index order which was provided by '--samples'.
    /// Note that the length of --samples and --sample-names do not need to match. When such is the case,
    /// default values are used for the missing names.
    /// 
    /// Example: '--sample 0-2 5 6 --sample-names ART04 ART16 ART20'{n}
    ///   - index  0      1      2      5     6{n}
    ///   - name   ART04  ART16  ART20  Ind5  Ind6{n}
    #[clap(short='n', long, multiple_values(true))]
    pub sample_names: Vec<String>,

    /// Output directory where results will be written.
    /// 
    /// Note that grups-rs will create the specified leaf directory if it is not present, by does not 
    /// allow itself from creating parent directories.
    #[clap(short, long, default_value("grups-output"), parse(try_from_os_str=valid_output_dir))]
    pub output_dir: PathBuf,

    /// Overwrite existing output files.
    /// 
    /// By default, grups-rs does not allow itself from overwriting existing results files. Use this flag
    /// to force this behaviour.
    #[clap(short='w', long)]
    pub overwrite: bool,
}

/// Estimate the raw average genetic PairWise Differences between individuals
/// 
/// Estimate the observed average pairwise mismatch rate (or 'PairWise Differences') between individuals 
/// within a pileup file.
#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug, Default, Serialize, Deserialize)]
pub struct PwdFromStdin {
    /// Enable self-comparison mode on individuals.
    /// 
    /// Note that a minimal depth of 2 is required when comparing individuals to themselves, since at least two
    /// alleles are required to effectively compute pairwise differences at a given site. This can be activated
    /// by specifying a library-wise value of --min-depth>=2 for all samples.
    /// 
    /// Note that in cases where the specified --min-depth was set to 1 for a given individual, grups-rs will
    /// automatically temporarily rescale the given value of --min-depth to 2 when performing self-comparisons.
    /// 
    #[clap(short='S', long)]
    pub self_comparison: bool,
    /// Filter out tri-allelic sites when given a list of SNP-coordinates.
    /// 
    /// Note that this arguement requires the use of --targets to provide the program with a list of known 
    /// coordinates. The given coordinate file must furthermore explicitely provide with REF/ALT alleles.
    /// 
    /// See the --targets argument for additional information regarding valid formats.
    #[clap(short='k', long)]
    pub known_variants: bool,

    /// Do not perform comparison, but rather print out the pileup lines where a comparison could have been made.
    ///
    /// This is a legacy argument from the initial implementation of grups. Note that lines are printed as soon
    /// as there is a valid comparison between ANY pair of individuals. Thus, applying sites filtration may not
    /// prove very effective when using --self-comparison.
    #[clap(short='f', long)]
    pub filter_sites: bool,

    /// Consider deletions when performing comparisons.
    /// 
    /// By default, deletion and missing characters ('*') are not counted as selectible nucleotide for comparison.
    /// Using this flag will instead mark them as valid matching positions.
    #[clap(long, action(ArgAction::SetFalse))]
    pub consider_dels: bool,

    /// Do not print jackknife blocks for each individual
    /// 
    /// By default, grups-rs will keep track of the pairwise mismatch rate within windows of size '--blocksize'.
    /// This information can prove useful to investigate PMR values in sliding windows (this can be visualized
    /// with the 'grups.plots' companion shiny interface), or to compute jackknife estimates of variance.
    /// 
    /// Using --no-print-block will prevent grups-rs from computing average pairiwise differences in non-overlapping
    /// blocks
    /// 
    #[clap(short='J', long)]
    pub no_print_blocks: bool,

    /// Change the window size of each jackknife block.
    /// 
    /// Note that lower block values may drastically increase the memory footprint of grups-rs.
    #[clap(short='b', long, required(false), default_value("1000000"))]
    pub blocksize: u32,
    
    /// Provide with the minimal sequencing depth required to perform comparison.
    /// 
    /// Note that the length of the provided vector does not need to be the same as the one
    /// provided for individuals. When such is the case, values of --min-depth will "wrap" 
    /// around, and start again from the beginning for the next individuals.
    /// 
    /// Example:  '--samples 0-4 --min-depth 2 3 5 ' will result in:{n}
    ///   - Ind  : 0 1 2 3 4{n}
    ///   - Depth: 2 3 5 2 3{n}
    #[clap(short='x', long, multiple_values(true), required(false), default_values(&["1","1"]))]
    pub min_depth: Vec <u16>,

    /// 0-based column index of the individuals that should be compared
    /// 
    /// Argument may accept slices (inclusive), such as --samples 0-3 and/or discrete integer values
    /// such as --samples 7 8. 
    /// 
    /// Example:{n}
    /// '--samples 0-3 7 8' will be parsed as: [0, 1, 2, 3, 7, 8]{n}
    /// 
    /// Note that samples are specified by 0-based index value. Thus, use --samples 0-3 if you have
    /// 4 samples within your pileup file, and wish to investigate all pairwise comparisons between 
    /// them.
    #[clap(short='s', long, multiple_values(true), default_values(&["0", "1"]))]
    pub samples: Vec<String>,

    /// Exclude transitions from the input targets file.
    /// 
    /// Note that this arguement requires the use of --targets to provide the program with a list of coordinates.
    /// The given coordinate file must furthermore explicitely provide with the REF/ALT alleles. See the --targets
    /// argument for additional information regarding valid file formats.
    #[clap(long)]
    pub exclude_transitions: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Serialize, Deserialize)]
pub enum Mode {
    Vcf,
    Fst,
    FstMmap
}

impl Default for Mode {
    fn default() -> Self {Self::Vcf}
}

#[derive(Debug, Copy, Clone, ArgEnum, Serialize, Deserialize)]
pub enum RelAssignMethod {Zscore, SVM }

impl Default for RelAssignMethod {
    fn default() -> Self {Self::SVM}
}

impl Display for RelAssignMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zscore => write!(f, "Z-score assignation"),
            Self::SVM    => write!(f, "Ordinally Partitionned SVM fitting")
        }
    }
}

/// Estimate genetic relatedness through pedigree simulations.
/// 
/// Perform genetic relatedness estimation between pileup individuals, by first running the pwd-from-stdin 
/// module, and performing pedigree simulations in one go.
#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug, Default, Serialize, Deserialize)]
pub struct PedigreeSims {
    /// Proportion of SNPs to keep at true frequency (in percentage).
    /// 
    /// grups-rs will use the specified --af-downsampling-rate as a probability to randomly
    /// simulate allele fixation during pedigree simulations.
    /// 
    /// Note that rates are specified as percentages.
    /// 
    #[clap(short='d', long, default_value("0.0"), parse(try_from_str=percent_str_to_ratio))]
    pub af_downsampling_rate: f64,

    /// Proportion of filtered SNP positions to include in the analysis (in percentage).
    /// 
    /// grups-rs will use the specified --snp-downsampling-rate as a probability to randomly
    /// ignore positions during pedigree simulations.
    /// 
    /// Note that rates are specified as percentages.
    #[clap(short='D', long, default_value("0.0"), parse(try_from_str=percent_str_to_ratio))]
    pub snp_downsampling_rate: f64,

    /// Contamination rates (or rate ranges) for each pileup individual. 
    /// 
    /// The provided argument(s) may accept hard set values, such as '--contam-rate 3.0', or ranges, such as '--contam-rate 0.0-5.0'.
    /// When a contamination rate is provided in the form of a range, pedigree-specific values are picked from a uniform distribution,
    /// within that defined range.
    /// 
    /// Note that contamination rates are specified as percentage, and are tied to the individuals contained within the input pileup.
    /// Specifying a constant rate for all individuals is also possible, by defining a unique value or range.
    /// 
    /// Example:{n}
    /// 
    /// grups-rs pedigree-sims [...] --samples 0-5 --contam-rate 2.0{n}
    /// 
    /// ...implies that all five pileup indidivuals will be assigned a set contamination rate of 2%, during pedigree simulations.{n}
    /// 
    /// In general, keep in mind that contamination rate values are recycled, if the number of specified values is lower than the number of
    /// examined pileup individuals.
    /// 
    #[clap(short='Q', long, required(false), multiple_values(true), default_values(&["0", "0"]), parse(try_from_str=parse_pedigree_param))]
    pub contamination_rate: Vec<Vec<f64>>,

    /// Sequencing error rates (or rate ranges) for each pileup individual.
    /// 
    /// The provided argument(s) may accept hard set values, such as '--seq-error-rate 1.0', or ranges, such as '--seq-error-rate 1.0-3.0'.
    /// When a sequencing error rate is provided in the form of a range, pedigree-specific values are picked from a uniform distribution,
    /// within that defined range.
    /// 
    /// Note that sequencing error rates are specified as percentage, and are tied to the individuals contained within the input pileup.
    /// Specifying a constant rate for all individuals is also possible, by defining a unique value or range.
    /// 
    /// Example:{n}
    /// 
    /// grups-rs pedigree-sims [...] --samples 0-7 --seq-error-rate 1.0{n}
    /// 
    /// ...implies that all seven pileup indidivuals will be assigned a set sequecing error rate of 1%, during pedigree simulations.{n}
    /// 
    /// In general, keep in mind that sequencing error rate values are recycled if the number of specified values is lower than the number of
    /// examined pileup individuals.
    #[clap(short='U', long, required(false), multiple_values(true), parse(try_from_str=parse_pedigree_param))]
    pub seq_error_rate : Option<Vec<Vec<f64>>>,

    /// Path to a directory containing a database of phased modern human genotypes, such as the 1000g-phase3 dataset.
    /// 
    /// This database may be in the form of a set of VCF files (default). In that case, grups-rs will look for, and 
    /// load any file ending with the .vcf[.gz] file extension.
    /// 
    /// Alternatively, users may use a set of FSA-encoded files. To use FSA-encoded files, users must explicitly specify
    /// the use of this type of data input, with the --mode argument. When specified, grups-rs will search for any file
    /// ending with the .fst[.frq] file extension. 
    /// 
    /// See the documentation of the 'fst' module, for a more detailled explanation regarding how to generate FSA-encoded
    /// databases.
    /// 
    #[clap(short='F', long, parse(try_from_os_str=valid_input_directory))] // default_value(r#"./data/founders/1000G-phase3-v5a"#),
    pub data_dir: PathBuf,

    /// Define the expected data input type for pedigree simulations.
    /// 
    /// This argument is closely tied to the --data-dir, and will define which type of files should grups-rs look for, as well
    /// as how to load them into memory. 
    /// 
    /// When using '--mode fst-mmap' grups-rs will search for files ending with the .fst[.frq] file extention
    /// 
    /// When using '--mode fst', grups-rs will search for files ending with the .fst[.frq] file extention. These files are then
    /// sequentially loaded into RAM and queried. The 'fst' mode may be faster if you have a lot of comparisons to apply, but has
    /// has a higher memory footprint compared to the 'fst-mmap' mode. Thus, it is only recommended if your input .fst[.frq] files
    /// are located on an HDD drive.
    /// 
    /// when using '--mode vcf' (default), grups-rs will search for files ending with the .vcf[.gz] extention, and will directly
    /// query these files. This mode is generally slower, but may yield higher performance if you have highly covered data, and
    /// your vcf files were previously pre-filtered to only contain bi-allelic SNPs.
    /// 
    #[clap(short='I', long, arg_enum, default_value("vcf"))]
    pub mode: Mode,


    /// Path to a directory containing a set of chromosome-specific genetic recombination maps, such as the HapMap-phaseII dataset.
    /// 
    /// Note that GRUPS-rs will search for, and attempt to load any file ending with the '.txt' file extension within that directory.
    /// It is therefore highly recommended that this directory ONLY contains the required recombination map files, and nothing else.
    /// 
    /// The expected input is a headed, tab separated file with columns '<chr> <pos(bp)> <rate(cM/Mb)> <Map(cM)>'{n}
    /// Example:
    /// Chromosome  Position(bp)    Rate(cM/Mb)     Map(cM)
    /// chr22       16051347        8.096992        0.000000
    /// chr22       16052618        8.131520        0.010291
    /// ...         ...             ...             ...
    /// 
    /// 
    #[clap(short='G', long, required(false), parse(try_from_os_str=valid_input_directory))] // default_value("./data/recombination/GRCh37"),
    pub recomb_dir: PathBuf,

    /// Number of pedigree simulation replicates to perform for each pairwise comparisons.
    /// 
    /// The default provided value of 100 should be considered a bare-minimum. Value in the range 500 to 1000 replicates is 
    /// recommended.
    #[clap(short='R', long, default_value("100"))]
    pub reps: u32,

    /// Minimum required (super-)population allele frequency to include a given SNP during simulations.
    /// 
    /// Note that grups-rs will re-compute the observed average PWD after simulations, to excluding and correct for
    /// any position previously excluded through this threshold.
    #[clap(short='m', long, default_value("0.00"))]
    pub maf: f32,

    /// Source population used during pedigree simulations.
    /// 
    /// Source (super-)population from which founder individuals are selected during pedigree simulations.
    /// 
    /// Note that when using '--mode vcf', grups-rs may only use populations for which a <POP>_AF annotation is present
    /// in each INFO field of the VCF files used. To generate finer grained population allele frequencies, we recommend
    /// either the use of the 'bcftools +fill-tags' plugin, or the use of the '--compute-pop-afs' argument when generating
    /// FSA-encoded dataset with the 'grups-rs fst' module.
    #[clap(short='P', long, default_value("EUR"))]
    pub pedigree_pop: String,

    /// Contaminating population used during pedigree simulations.
    /// 
    /// (Super-)population from which contaminating individuals are selected during pedigree simulations. 
    /// 
    /// Note that when using '--mode vcf', grups-rs may only use populations for which a <POP>_AF annotation is present
    /// in each INFO field of the VCF files used. To generate finer grained population allele frequencies, we recommend
    /// either the use of the 'bcftools +fill-tags' plugin, or the use of the '--compute-pop-afs' argument when generating
    /// FSA-encoded dataset with the 'grups-rs fst' module.
    #[clap(short='C', long, multiple_values(true), default_values(&["EUR"]))]
    pub contam_pop: Vec<String>,

    /// Number of random individual genomes with which to contaminate pedigree simulations.
    /// 
    /// Format: <INT> <INT> ...
    /// Default: Contaminate all pedigree individuals with a single contaminating individual.
    #[clap(short='N', long, multiple_values(true), default_values(&["1"]))]
    pub contam_num_ind: Vec<usize>,

    /// Path to input pedigree definition file.
    #[clap(short='T', long, required(false), parse(try_from_os_str=valid_input_file))] // default_value(r#"./data/pedigrees/default.ped"#),
    pub pedigree: PathBuf,
    
    // Path to input Panel Definition files.
    #[clap(short='p', long, parse(try_from_os_str=valid_input_file))]
    pub panel: Option<PathBuf>,

    /// Number of additional parallel decompression threads.
    /// 
    /// Can increase performance when working with BGZF compressed .vcf.gz files. Note that this parameter has no effect when working with
    /// uncompressed .vcf or and .fst[.frq] files.
    /// 
    #[clap(short='#', long, default_value("0"))]
    pub decompression_threads: usize,

    /// Provide the RNG with a set seed.
    #[clap(long, required(false), default_value_t=fastrand::u64(u64::MIN..=u64::MAX))]
    pub seed: u64,

    /// Select the method for most likely relationship assignment
    /// 
    /// zscore: Perform minimum zscore assignation, i.e. the distribution with the lowest z-score from the observed PWD is selected as the most likely candidate.
    /// Computationally inexpensive, but can provide with surprising results, when the different distributions carry drastically different standard deviations.
    /// 
    /// svm: Compute treshold using Ordinally Partitionned Support Vector Machines. Binary SVMs are instantiated sequentially, from the lowest relatedness order 
    /// to the highest, in terms of average PWD. Each SVM is fitted against the hypothesis that the observed PWD belongs to a higher degree than given distribution.
    /// i.e., the question asked by the svm can be translated into: "Is this observed PWD greater than than the currently observed distribution?"
    /// The most likely relationship is assigned as soon as an SVM answers 'no'.
    #[clap(long, arg_enum, default_value("svm"))]
    pub assign_method: RelAssignMethod
}

/// Convert VCF files into sets of FSA-encoded indexes.
/// 
/// Encode large sets of phased VCF reference files in the form of Deterministic Acyclic Finite State Automata 
/// to enable a higher query performance when performing pedigree simulations.
#[derive(Args, Debug, Default, Serialize, Deserialize)]
pub struct VCFFst {
    /// Path to input Panel Definition files.
    #[clap(short='p', long, parse(try_from_os_str=valid_input_file))]
    pub panel: Option<PathBuf>,

    /// Population Subset
    /// 
    /// Subset the index by a given number of (super)-population. (e.g. EUR, AFR). Note that at least one pedigree population and one contaminating population
    /// are required to obtain valid index files.
    #[clap(short='P', long, multiple_values(true))]
    pub pop_subset: Option<Vec<String>>,

    /// Output directory
    #[clap(short='o', long, parse(try_from_os_str=valid_output_dir))]
    pub output_dir: PathBuf,

    /// Path to input VCF genomes for founder individuals (1000G)
    #[clap(short='d', long, parse(try_from_os_str=valid_input_directory))] // default_value(r#"./data/founders/1000G-phase3-v5a"#),
    pub vcf_dir: PathBuf,

    /// Number of parallel CPU processes when performing FST-Indexation
    /// 
    /// Parallelization is dispatched according to the number of separate vcf(.gz) files.
    #[clap(short='@', long, default_value("1"))]
    pub threads: usize,

    ///Number of additional parallel decompression threads when decompressing (BGZF compressed files only).
    /// 
    /// Parallelization is dispatched according to the number of replicates.
    #[clap(short='#', long, default_value("0"))]
    pub decompression_threads: usize,

    /// Recalculate allele frequencies for each (super-)population
    /// 
    /// Recompute population allele frequencies for each population and super-population tag that can be found within the input
    /// panel definition file. 
    /// 
    /// For some publicly available datasets, such as the 1000g-phase3 adding this flag can allow for the use of smaller populations as reference
    /// or to use customly defined populations.
    /// 
    /// When unspecified, the program will instead look for <POP>_AF tags within the VCF's INFO field. These tags can be generated
    /// Using the bcftools '+fill-tags' plugin.
    #[clap(short='F', long)]
    pub compute_pop_afs: bool

}

impl PwdFromStdin {
    /// Sanity check : depth must indeed be > 2 when performing self-comparison.
    ///
    /// # Errors
    ///  if any of the values provided through `--min-depth` is lower than 2 while `--self-comparison` is on.
    /// 
    /// # @TODO:
    /// - put this in `Comparison::new()` ?? -> This would allow mix-matching batch mode and self-comparison,
    ///   but could be a bit confusing for users..
    pub fn check_depth(&self) -> Result<(), ParserError> {
        if self.min_depth.iter().any(|&x| x < 1) {        
            return Err(ParserError::InsufficientDepthError)
        }
        Ok(())
    }
}

/// Command line interface argument parser.
/// 
/// TODO : - `get_results_file_prefix()` and `get_blocks_output_files()` should not be the responsability
///          of this struct. --> migrate to `pwd_from_stdin::io.rs`
///        - add deserialization method. Users could thus fully reproduce a previous run with ease. keep it FAIR. 
impl Common {
    /// Sanity Check: The program should leave if the user did not provide any pileup input, either through
    /// `--pileup` or through stdinput. Without this, our program would wait indefinitely for the stdin buffer.
    /// 
    /// # Errors
    /// - if the user did not provide an input file, neither from stdin, nor through the `--pileup` argument.
    pub fn check_input(&self) -> Result<(), ParserError> {
        if atty::is(atty::Stream::Stdin) && self.pileup.is_none() {
            return Err(ParserError::MissingPileupInput)
        }
        Ok(())
    }

    /// Get a generic filename for our output files. If the user used `--pileup`, this will become its file stem.
    /// If the user used stdin, this will become a generic name -> "pwd_from_stdin-output" 
    ///
    /// # @TODO: This function should be the one responsible of defining the default filename. Stay dry.
    /// 
    /// # Errors 
    /// - if a default file-prefix cannot be created from the input pileup filestem.
    pub fn get_file_prefix(&self, subdir: Option<&str>) -> Result<PathBuf> {
        let default_prefix = String::from("pwd_from_stdin-output");
        let file_prefix = Path::new(self.pileup.as_ref()
            .unwrap_or(&default_prefix))
            .file_stem()
            .ok_or_else(||anyhow!(ParserError::ParseOutputPrefix))
            .loc( "While parsing command line arguments" )?;

        // Get the path/filename-prefix of all of our outputs.
        // @TODO: Don't push things like that, it's yucky. 
        let mut parsed_file = PathBuf::new();
        parsed_file.push(&self.output_dir);
        parsed_file.push(subdir.unwrap_or(""));
        parsed_file.push(file_prefix);
        Ok(parsed_file)
    }

    /// Check if a given file already exists ; raise an error if such is the case, and the user did not explicitly 
    /// allow file overwriting.
    /// 
    /// # Errors
    /// - If the provided `pathbuf` already exists and the user did not specifically allow for file
    ///   overwrite using the `--overwrite` argument
    /// 
    /// # Panics
    /// - if the provided `pathbuf` fails to get parsed as a string.
    pub fn can_write_file(&self, pathbuf: &Path) -> Result<bool> {
        if ! self.overwrite && pathbuf.exists() {   // Check if this file already exists and/or if overwrite is allowed.
            return Err(ParserError::CannotOverwrite(pathbuf.display().to_string()))
                .loc( "While parsing command line arguments" )
        }
        Ok(true)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum FileEntity {File, Directory}

impl Display for FileEntity {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            Self::File      => write!(f, "File"),
            Self::Directory => write!(f, "Directory"),
        }
    }
}

impl FileEntity {
    fn validate(&self, path: &Path) -> Result<(), ParserError> {
        use ParserError::InvalidFileEntity;
        let valid = match self {
            Self::File      => path.is_file(),
            Self::Directory => path.is_dir()
        };

        if valid {
            Ok(())
        } else {
            Err(InvalidFileEntity(*self, path.display().to_string()))
        }
    }
}

fn assert_filesystem_entity_is_valid(s: &OsStr, entity: &FileEntity) -> Result<()> {
    use ParserError::MissingFileEntity;
    let path = Path::new(s);
    if ! path.exists() {
        return Err(MissingFileEntity(*entity, path.display().to_string()))
            .loc("While parsing arguments.")
    }

    entity.validate(path).loc("While parsing arguments.")
}

fn valid_input_directory(s: &OsStr) -> Result<PathBuf> {
    assert_filesystem_entity_is_valid(s, &FileEntity::Directory)
        .loc("While checking for directory validity")?;
    Ok(PathBuf::from(s))
}

fn valid_input_file(s: &OsStr) -> Result<PathBuf> {
    assert_filesystem_entity_is_valid(s, &FileEntity::File)
        .loc("While checking for file validity")?;
    Ok(PathBuf::from(s))
}

fn valid_output_dir(s: &OsStr) -> Result<PathBuf> {
    if ! Path::new(s).exists() {
        std::fs::create_dir(s)?;
    }
    assert_filesystem_entity_is_valid(s, &FileEntity::Directory)
        .loc("While checking for directory validity")?;
    Ok(PathBuf::from(s))
}

/// Convert a user-defined string "range" into a vector of integers.
/// "9-14" thus becomes [9, 10, 11, 12, 13, 14]
/// Note that the range is fully inclusive. 
fn parse_user_range<T>(s: &str) -> Result<Vec<T>, <T as FromStr>::Err> 
where   T       : FromStr + Add<Output = T> + Ord + One,
        Range<T>: Iterator<Item = T>,
{
    match s.split_once('-') {
            Some(t) => Ok((t.0.parse::<T>()?..t.1.parse::<T>()?+One::one()).collect::<Vec<T>>()),
            None    => Ok(vec![s.parse::<T>()?])
    }
}

fn percent_str_to_ratio(s: &str) -> Result<f64>{
    use ParserError::ParseRatio;

    const MIN_PERCENT: f64 = 0.0;
    const MAX_PERCENT: f64 = 100.0;

    let percent = s.parse::<f64>()?;

    // Ensure the user input lies between the [0% - 100%] range.
    match percent < MAX_PERCENT || percent > MIN_PERCENT {
        true  => Ok(percent/100.0),
        false => Err(anyhow!(ParseRatio(MIN_PERCENT, MAX_PERCENT))).with_loc(|| format!("While parsing {s}"))
    }
}

fn parse_pedigree_param(s: &str) -> Result<Vec<f64>> {
    use ParserError::ParseRange;
    let vec: Vec<f64> = s.split('-')
        .map(percent_str_to_ratio)
        .collect::<Result<Vec<f64>, _>>()
        .with_loc(|| format!("While parsing the provided string: {s}"))?;

    if vec.is_empty() || vec.len() > 2 {
        return Err(anyhow!(ParseRange("Found multiple dashes")))
        .with_loc(|| format!("While parsing the provided string: {s}"))
    }

    Ok(vec)
}

/// Convert a vector of Strings with user-input ranges to a single vector of integers.
///  - Input will most likely stem from the command line parser, where users are not expected
///    to write every single value they would like to input.
/// 
///     --> ["1-6", "8"] for the user, becomes [1, 2, 3, 4, 5, 6, 8] for our program.
/// 
///  - Return a vector of generic integers, boxed within a Result.
///
/// # @TODO:
///   - Please make this more elegant? We could use a `flat_map` at some point...
/// 
/// # Example
///```
/////use pwd_from_stdin::parser::parse_user_ranges;
/////let user_input:   Vec<&str> = vec!["5", "1-3", "7"];
/////let parsed_input: Vec<usize>  = parse_user_ranges(&user_input, "chr").expect("error");
/////assert_eq!(parsed_input, vec![1, 2, 3, 5, 7])
///```
/// 
/// # Errors
///  returns a `ParseError` if the provided ranges cannot be parsed into integers.
pub fn parse_user_ranges<T>(ranges: &[String], arg: &str) -> Result<Vec<T>, ParserError>
where   T                   : FromStr + Add<Output = T> + Ord + One,
        Range<T>            : Iterator<Item = T>,
        <T as FromStr>::Err : ToString,

{
    let parsed_ranges = ranges.iter()
        .map(|s| parse_user_range(s))
        .collect::<Result<Vec<Vec<T>>, _>>();

    // @TODO - This is inefficient: We're  performing a two-pass iteration while we could 
    //         flatten the structure immediately. use .flat_map() 
    let mut parsed_ranges: Vec<T> = match parsed_ranges {
        Ok(vec) => vec.into_iter().flatten().collect(),
        Err(err) => return Err(ParserError::ParseArg{arg: arg.to_string(), err: err.to_string()}),
    };
    parsed_ranges.sort();
    parsed_ranges.dedup();
    Ok(parsed_ranges)
}