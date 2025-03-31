extern crate parser;
extern crate logger;

use std::fs::File;

use parser::{Cli, Commands::*};
use genome::Genome;

#[macro_use]
extern crate log;

use anyhow::{anyhow, Result};

pub fn cite() {
    // If this ever becomes bloated, consider using the 'indoc' crate.
    const CITATIONS: &str = r###"
    A. If you plan to use GRUPS-rs in your work, please cite the corresponding 
       paper of GRUPS-rs (1), as well as the original publication of GRUPS (2):

        1. Lefeuvre M, Martin M, Jay F, Marsolier M, Bon C. GRUPS-rs, a
        high-performance ancient DNA genetic relatedness estimation software
        relying on pedigree simulations. Hum Popul Genet Genom 2024;4(1):0001.
        https://doi.org/10.47248/hpgg2404010001 

        2. Martin MD, Jay F, Castellano S, Slatkin M. Determination of
        genetic relatedness from low-coverage human genome sequences using
        pedigree simulations. Mol Ecol. 2017; 26: 4145– 4157.
        https://doi.org/10.1111/mec.14188 


    B. If you plan to use the HapMap-phaseII-b37 recombination maps described
       in the README and examples, please cite the International HapMap
       consortium paper(s):

       1. The International HapMap Consortium. A second generation human 
          haplotype map of over 3.1 million SNPs. Nature 449, 851–861 (2007).
          https://doi.org/10.1038/nature06258

       2. Kong, A., Gudbjartsson, D., Sainz, J. et al. A high-resolution 
          recombination map of the human genome. Nat Genet 31, 241–247 (2002).
          https://doi.org/10.1038/ng917
    

    C. If you plan to use the 1000genomes-phase 3 variants callset described in
       the README and examples, please cite the 1000 Genomes Project Consortium

       1. The 1000 Genomes Project Consortium. A map of human genome variation
          from population-scale sequencing. Nature 467, 1061–1073 (2010).
          https://doi.org/10.1038/nature09534

    "###;
    println!("{CITATIONS}");
}

pub fn run(cli: Cli) -> Result<()> {
    // ----------------------------- Set seed (randomly assigned by parser-rs if none was provided.)
    if let PedigreeSims{ref ped, common: _, pwd:_ } = cli.commands {
        fastrand::seed(ped.seed);
        trace!("Seed: {}", fastrand::get_seed());
    }
    // ----------------------------- Initialize genome.
    let genome = match cli.commands {
        PedigreeSims { ref common, pwd: _, ped: _} | PwdFromStdin { ref common, pwd: _ } => {
            info!("Indexing reference genome...");
            let mut genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None       => Genome::default(),
            };
            // ---- Split genome into autosome and X-chromosome chunks.            
            let x_chromosome = genome.pop_xchr();

            if common.x_chromosome_mode {
                genome = x_chromosome.ok_or_else(||anyhow!("Missing Xchromosome"))?;
            };
            Some(genome)
        },
        _ => None
    };

    // ---- Run with autosomes.
    _run(&cli, genome)?;

    Ok(())
}
pub fn _run(cli: &Cli, genome: Option<Genome>) -> Result<()> {
    match &cli.commands {
        PedigreeSims {common, pwd, ped} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            let genome = genome.ok_or_else(|| anyhow!("Error: Missing genome"))?; // @TODO better handling

            // ----------------------------- Run PWD_from_stdin.
            let (mut comparisons, output_files) = pwd_from_stdin::run(common, pwd, &requested_samples, &genome)?;
            comparisons.write_pwd_results(!pwd.no_print_blocks, &output_files)?;

            // ----------------------------- Run Pedigree-sims
            pedigree_sims::run(common, ped, &requested_samples, &mut comparisons)?;

        },
        
        PwdFromStdin {common, pwd} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            let genome = genome.ok_or_else(|| anyhow!("Error: Missing genome"))?; // @TODO better handling

            // ----------------------------- Run PWD_from_stdin.
            let (comparisons, output_files) = pwd_from_stdin::run(common, pwd, &requested_samples, &genome)?;
            if ! pwd.filter_sites {
                comparisons.write_pwd_results(!pwd.no_print_blocks, &output_files)?;
            }
        },

        FST {fst: fst_cli} => {
            vcf_fst::run(fst_cli)?
        },

        FromYaml{yaml} => {
            let yaml_file = File::open(yaml)?;
            let cli: Cli = match serde_yaml::from_reader(yaml_file){
                Ok(cli)  => cli,
                Err(e)   => return Err(anyhow!("Unable to deserialize arguments from {yaml:?} file: [{e}]"))
            };
            self::run(cli)?;
        },

        Cite => {
            cite();
        }
    };
    Ok(())
}