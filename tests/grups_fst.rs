use clap::Parser;

#[cfg(test)]
mod common;
use common::Fixture;
use common::test_grups_run;

use anyhow::Result;

#[test]
fn test_grups_fst() -> Result<()> {
    let vcf_dir    =  Fixture::copy("vcf/binary-2FIN-1ACB-virtual/");
    let output_dir  = Fixture::blank("test_grups_fst");

    let args = format!("grups fst
        --vcf-dir    {vcf_dir}
        --output-dir {output_dir}
    ");

    eprintln!("{args}");

    let cli = parser::Cli::parse_from(args.split_whitespace());
    grups::run(cli)?;
    
    // Ensure grups can correctly estimate results using the newly generated FST index.
    test_grups_run(parser::Mode::Fst, &format!("{output_dir}"));

    Ok(())

}