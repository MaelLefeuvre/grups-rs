extern crate parser;
extern crate logger;

use parser::{Cli, Commands::*};
use genome::Genome;

#[macro_use]
extern crate log;

use std::error::Error;

pub fn cite() {
    // If this ever becomes bloated, consider using the 'indoc' crate.
    const CITATIONS: &str = r###"
    A. If you plan to use GRUPS-rs in your work, please cite the original
       publication of GRUPS:

        1. Martin, MD, Jay, F, Castellano, S, Slatkin, M. Determination of
        genetic relatedness from low-coverage human genome sequences using
        pedigree simulations. Mol Ecol. 2017; 26: 4145– 4157.
        https://doi.org/10.1111/mec.14188 


    B. If you plan to use the HapMap-II-b37 recombination maps described in the
       README and examples, please cite the International HapMap consortium:

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

/// @TODO : Stay dry...
pub fn run(cli: Cli) -> Result<(), Box<dyn Error>> {
    match cli.commands {
        PedigreeSims {common, pwd, ped} => {
            // ----------------------------- Set seed (randomly assigned by parser-rs if none was provided.)
            fastrand::seed(ped.seed);
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None => Genome::default(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (mut comparisons, output_files) = pwd_from_stdin::run(&common, &pwd, &requested_samples, &genome)?;
            comparisons.write_pwd_results(pwd.print_blocks, &output_files)?;
            pedigree_sims::run(common, *ped, &mut comparisons)?;

        },
        
        PwdFromStdin {common, pwd} => {
            // ----------------------------- Parse Requested_samples
            let requested_samples: Vec<usize> = parser::parse_user_ranges(&pwd.samples, "samples")?;
            // ----------------------------- Initialize genome.
            info!("Indexing reference genome...");
            let genome = match &common.genome {
                Some(file) => Genome::from_fasta_index(file)?,
                None => Genome::default(),
            };
            // ----------------------------- Run PWD_from_stdin.
            let (comparisons, output_files) = pwd_from_stdin::run(&common, &pwd, &requested_samples, &genome)?;
            if ! pwd.filter_sites {
                comparisons.write_pwd_results(pwd.print_blocks, &output_files)?;
            }
        },

        FST {fst: fst_cli} => {
            vcf_fst::run(&fst_cli)?
        },

        FromYaml{yaml} => {
            let cli: Cli = match serde_yaml::from_reader(std::fs::File::open(&yaml)?) {
                Ok(cli)  => cli,
                Err(e) => return Err(format!("Unable to deserialize arguments from {yaml:?} file: [{e}]").into())
            };
            self::run(cli)?;
        },

        Cite => {
            cite();
        }
    };
    Ok(())
}