use clap::Parser;

#[cfg(test)]
mod common;
use common::Fixture;
use common::test_grups_run;

#[test]
fn test_grups_fst() {
    let vcf_dir    =  Fixture::copy(&format!("vcf/binary-2FIN-1ACB-virtual/"));
    let output_dir  = Fixture::blank(&format!("test_grups_fst"));

    let args = format!("grups fst
        --vcf-dir    {vcf_dir}
        --output-dir {output_dir}
    ");

    let cli = parser::Cli::parse_from(args.split_whitespace());
    grups::run(cli).unwrap();
    
    // Ensure grups can correctly estimate results using the newly generated FST index.
    test_grups_run(parser::Mode::Fst, &format!("{output_dir}"));

}