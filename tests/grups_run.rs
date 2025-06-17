mod common;
use common::GrupsRunnerBuilder;
#[cfg(test)] use pretty_assertions::assert_eq;

#[test]
fn test_grups_run_vcf_autosomes() {
    let runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir("vcf/binary-2FIN-1ACB-virtual-autosomes/")
        .set_pileup("pileup/parents-offspring.pileup")
        .set_recomb_dir("recombination-maps/autosomes")
        .set_pedigree("pedigree/tiny_pedigree.txt")
        .set_output_dir("grups-test-output")
        .set_contam_pop("AFR")
        .set_samples("0-2")
        .set_mode(parser::Mode::Vcf)
        .overwrite()
        .build()
        .unwrap();

    runner.run();
    validate_file!("expect/parents-offspring.pwd", runner.output_pwd().unwrap());
    validate_file!("expect/parents-offspring-vcf.result", runner.output_results().unwrap());
    //validate_file!("expect/parents-offspring-vcf.probs", runner.output_probs().unwrap());
}

#[test]
fn test_grups_run_fst_autosomes() {
    let runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir("fst/binary-2FIN-1ACB-virtual-autosomes/")
        .set_pileup("pileup/parents-offspring.pileup")
        .set_recomb_dir("recombination-maps/autosomes")
        .set_pedigree("pedigree/tiny_pedigree.txt")
        .set_output_dir("grups-test-output")
        .set_contam_pop("AFR")
        .set_samples("0-2")
        .set_mode(parser::Mode::Fst)
        .overwrite()
        .build()
        .unwrap();

    runner.run();
    validate_file!("expect/parents-offspring.pwd", runner.output_pwd().unwrap());
    validate_file!("expect/parents-offspring-fst.result", runner.output_results().unwrap());
    //validate_file!("expect/parents-offspring-fst.probs", runner.output_probs().unwrap());
}

#[test]
fn test_grups_run_vcf_chr_x() {
    let runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir("vcf/binary-2FIN-1ACB-virtual-chrX/")
        .set_pileup("pileup/ash128-ash133.virtual.min-depth-5.chrX.pileup")
        .set_recomb_dir("recombination-maps/chrX")
        .set_pedigree("pedigree/chrX-pedigree.ped")
        .set_output_dir("grups-test-output")
        .set_targets("targets/chrX-maf10.subset.filtered.tsv")
        .set_contam_pop("AFR")
        .set_samples("0-1")
        .set_mode(parser::Mode::Vcf)
        .x_chromosome_mode()
        .exclude_transitions()
        .set_reps(100)
        .set_seed(42)
        .overwrite()
        .build()
        .unwrap();

    runner.run();

    validate_file!("expect/chrX-result.pwd", runner.output_pwd().unwrap());
    validate_file!("expect/chrX-result-vcf.result", runner.output_results().unwrap());
    //validate_file!("expect/chrX-result-vcf.probs", runner.output_probs().unwrap());
}

#[test]
fn test_grups_run_fst_chr_x() {
    let runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir("fst/binary-2FIN-1ACB-virtual-chrX/")
        .set_targets("targets/chrX-maf10.subset.filtered.tsv")
        .set_pileup("pileup/ash128-ash133.virtual.min-depth-5.chrX.pileup")
        .set_recomb_dir("recombination-maps/chrX")
        .set_pedigree("pedigree/chrX-pedigree.ped")
        .set_output_dir("grups-test-output-chrX")
        .set_contam_pop("AFR")
        .set_samples("0-1")
        .set_mode(parser::Mode::Fst)
        .set_reps(100)
        .set_seed(42)
        .x_chromosome_mode()
        .exclude_transitions()
        .overwrite()
        .build()
        .unwrap();

    runner.run();
    validate_file!("expect/chrX-result.pwd", runner.output_pwd().unwrap());
    validate_file!("expect/chrX-result-fst.result", runner.output_results().unwrap()); // Bugs out when running 'cargo llvm-cov --workspace' ...
    //validate_file!("expect/chrX-result-fst.probs", runner.output_probs().unwrap());
}

#[test]
fn test_grups_run_sex_specific_mode() {
    let runner = GrupsRunnerBuilder::new()
        .module("pedigree-sims")
        .set_data_dir("fst/binary-2FIN-1ACB-virtual-autosomes/")
        .set_pileup("pileup/parents-offspring.pileup")
        .set_recomb_dir("recombination-maps/autosomes")
        .set_pedigree("pedigree/chrX-pedigree.ped")
        .set_output_dir("grups-test-output")
        .set_contam_pop("AFR")
        .set_samples("0-1")
        .set_mode(parser::Mode::FstMmap)
        .sex_specific_mode()
        .overwrite()
        .build()
        .unwrap();

    runner.run();
    validate_file!(
        "expect/parents-offspring-sex-specific-sims/parents-offspring-Ind0-Ind1.sims",
        format!("{}-Ind0-Ind1.sims", runner.simulations_filestem())
    );
}