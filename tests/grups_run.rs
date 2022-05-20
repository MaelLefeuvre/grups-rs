mod common;

#[test]
fn test_grups_run_vcf() {
    let data_dir   = format!("vcf/binary-2FIN-1ACB-virtual/");
    common::test_grups_run(parser::Mode::Vcf, &data_dir);
}


#[test]
fn test_grups_run_fst() {
    let data_dir   = format!("fst/binary-2FIN-1ACB-virtual/");
    common::test_grups_run(parser::Mode::Fst, &data_dir);
}