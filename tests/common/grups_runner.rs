use clap::Parser;
use std::io::{BufRead, BufReader};
use std::fs::File;

use super::Fixture;

fn get_bufreader(filename: &str) -> BufReader<File> {
    let inner = File::open(filename)
        .expect(&format!("Failed to open test output file: {filename}"));
    BufReader::new(inner)
}

fn parse_line<'a, E: std::fmt::Debug>(line: &'a Result<String,E>, sep: &str) -> Vec<&'a str> {
    let line = line.as_ref().unwrap();
    line.split(sep).map(|field| field.trim()).collect::<Vec<&str>>()
}

fn assert_pwd_matches(
    filename: &str,
    expected_overlap: Vec<&str>,
    expected_avg_pwd: Vec<&str>,
){

    const OVERLAP_COL: usize = 1;
    const PWD_COL    : usize = 3;

    let result = get_bufreader(filename);

    for (i, line) in result.lines().skip(1).enumerate(){
        let line = parse_line(&line, " - ");
        assert_eq!(line[OVERLAP_COL], expected_overlap[i]);
        assert_eq!(line[PWD_COL], expected_avg_pwd[i]);
    }
}

fn assert_simulation_results(
    filename: &str,
    expected_relationships: Vec<&str>
) {

    const RELATIONSHIP_COL: usize = 1;

    let result = get_bufreader(filename);
    for (i, line) in result.lines().skip(1).enumerate(){
        let line = parse_line(&line, " - ");
        assert_eq!(line[RELATIONSHIP_COL], expected_relationships[i]);
        //assert_eq!(line[PWD_COL], expected_avg_pwd[i]);
    }
}

pub fn test_grups_run(mode: parser::Mode, data_dir: &str) {

    let mode_str = match mode {
        parser::Mode::Fst => "fst",
        parser::Mode::Vcf => "vcf",
    }; 

    const FILESTEM: &str = "parents-offspring";
    println!("before: {data_dir}");
    let data_dir   = Fixture::copy(&data_dir);
    println!("after: {data_dir}");

    //loop{}

    let pileup_dir = Fixture::copy(&format!("pileup/{FILESTEM}.pileup"));
    let recomb_dir = Fixture::copy("recombination-map/");
    let pedigree   = Fixture::copy("pedigree/tiny_pedigree.txt");
    let output_dir = Fixture::blank("grups-test-output");

    let args = format!("grups pedigree-sims
        --pileup {pileup_dir}
        --data-dir {data_dir} 
        --recomb-dir {recomb_dir}
        --pedigree {pedigree}
        --output-dir {output_dir}
        --overwrite
        --mode {mode_str}
        --samples 0-2
        --reps 10
    ",);

    println!("{args}");

    let cli = parser::Cli::parse_from(args.split_whitespace());
    grups::run(cli).unwrap();

    let output_pwd            = format!("{output_dir}/{FILESTEM}.pwd");
    let expected_overlap   = vec!["107"     , "86"      , "31"      ];
    let expected_avg_pwd   = vec!["0.18692" , "0.18605" , "0.29032" ];
    assert_pwd_matches(&output_pwd, expected_overlap, expected_avg_pwd);

    let output_results            = format!("{output_dir}/{FILESTEM}.result");
    let expected_relationships   = vec!["First Degree", "First Degree", "Unrelated"];
    assert_simulation_results(&output_results, expected_relationships);
}