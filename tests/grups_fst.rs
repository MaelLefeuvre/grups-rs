#[cfg(test)]
mod common;
use common::GrupsRunnerBuilder;

use anyhow::Result;
#[cfg(test)] use pretty_assertions::assert_eq;

#[test]
fn test_grups_fst() -> Result<()> {

    // ---- Create fst index
    let fst_runner = GrupsRunnerBuilder::new()
        .module("fst")
        .set_vcf_dir("vcf/binary-2FIN-1ACB-virtual-autosomes/")
        .set_output_dir("test_grups_fst")
        .build()
        .expect("GrupsRunner object should be buildable");

    fst_runner.run();
    // ---- Ensure grups can correctly estimate results using the newly generated FST index.
    let fst_dir = fst_runner.get_fixture("output-dir").to_string();
    let pedigree_sims_runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir(fst_dir.as_str())
        .set_pileup("pileup/parents-offspring.pileup")
        .set_recomb_dir("recombination-maps/autosomes/")
        .set_pedigree("pedigree/tiny_pedigree.txt")
        .set_output_dir("grups-test-output")
        .set_contam_pop("AFR")
        .set_samples("0-2")
        .set_mode(parser::Mode::Fst)
        .overwrite()
        .build()
        .expect("GrupsRunner Object should be buildable");
    pedigree_sims_runner.run();

    validate_file!("expect/parents-offspring.pwd", pedigree_sims_runner.output_pwd().unwrap());
    validate_file!("expect/parents-offspring-fst.result", pedigree_sims_runner.output_results().unwrap());
    Ok(())

}