use anyhow::{Result, bail};
use clap::Parser;
use std::collections::HashMap;
use std::path::Path;

use super::Fixture;

pub const TEST_PEDIGREE_SIMS_REPS: usize = 10;
pub const TEST_PEDIGREE_SIMS_SEED: usize = 42;

#[derive(Default)]
pub struct GrupsRunnerBuilder<'a> {
    module: Option<&'a str>,     // any
    pileup: Option<&'a str>,     // pwd-from-stdin + pedigree-sims
    targets: Option<&'a str>,    // pwd-from-stdin + pedigree-sims
    output_dir: Option<&'a str>, // pwd-from-stdin + pedigree-sims + fst
    samples: Option<&'a str>,    // pwd-from-stdin + pedigree-sims
    overwrite: bool,             // pwd-from-stdin + pedigree-sims
    x_chromosome_mode: bool,     // pwd-from-stdin + pedigree-sims
    exclude_transitions: bool,   // pwd-from-stdin + pedigree-sims
    pedigree: Option<&'a str>,   // pedigree-sims
    recomb_dir: Option<&'a str>, // pedigree-sims
    data_dir: Option<&'a str>,   // pedigree-sims
    mode: Option<parser::Mode>,  // pedigree-sims
    reps: Option<usize>,         // pedigree-sims
    contam_pop: Option<&'a str>, // pedigree-sims
    sex_specific_mode: bool,     // pedigree-sims
    seed: Option<usize>,         // pedigree-sims
    vcf_dir: Option<&'a str>,    // fst
    pop_subset: Option<&'a str>, // fst
    compute_pop_afs: bool,       // fst
}

impl <'a> GrupsRunnerBuilder<'a> {
    pub fn new() -> Self {
        Self{..Default::default()}
    }

    pub fn module(mut self, module: &'a str) -> Self {
        self.module = Some(module);
        self
    }

    pub fn set_mode(mut self, mode: parser::Mode) -> Self {
        self.mode = Some(mode);
        self
    }

    pub fn set_data_dir(mut self, path: &'a str) -> Self {
        self.data_dir = Some(path);
        self
    }

    pub fn set_recomb_dir(mut self, path: &'a str) -> Self {
        self.recomb_dir = Some(path);
        self
    }

    pub fn set_output_dir(mut self, path: &'a str) -> Self {
        self.output_dir = Some(path);
        self
    }

    pub fn set_pileup(mut self, path: &'a str) -> Self {
        self.pileup = Some(path);
        self
    }

    pub fn set_pedigree(mut self, path: &'a str) -> Self {
        self.pedigree = Some(path);
        self
    }

    #[allow(unused)]
    pub fn set_targets(mut self, path: &'a str) -> Self {
        self.targets = Some(path);
        self
    }

    pub fn set_samples(mut self, samples: &'a str) -> Self {
        self.samples = Some(samples);
        self
    }

    #[allow(unused)]
    pub fn set_reps(mut self, reps: usize) -> Self {
        self.reps = Some(reps);
        self
    }

    #[allow(unused)]
    pub fn set_seed(mut self, seed: usize) -> Self {
        self.seed = Some(seed);
        self
    }

    pub fn set_contam_pop(mut self, pop: &'a str) -> Self {
        self.contam_pop = Some(pop);
        self
    }

    pub fn overwrite(mut self) -> Self {
        self.overwrite = true;
        self
    }

    #[allow(unused)]
    pub fn x_chromosome_mode(mut self) -> Self {
        self.x_chromosome_mode = true;
        self
    }

    #[allow(unused)]
    pub fn exclude_transitions(mut self) -> Self {
        self.exclude_transitions = true;
        self
    }

    #[allow(unused)]
    pub fn sex_specific_mode(mut self) -> Self {
        self.sex_specific_mode = true;
        self
    }

    #[allow(unused)]
    pub fn set_vcf_dir(mut self, vcf_dir: &'a str) -> Self {
        self.vcf_dir = Some(vcf_dir);
        self
    }

    #[allow(unused)]
    pub fn set_pop_subset(mut self, pop_subset: &'a str) -> Self {
        self.pop_subset= Some(pop_subset);
        self
    }
    #[allow(unused)]
    pub fn compute_pop_afs(mut self) -> Self {
        self.compute_pop_afs = true;
        self
    }

    pub fn build(self) -> Result<GrupsRunner> {
        let mut args = vec!["grups-rs".to_string()];
        let mut filestem = None;
        let mut output_dir = String::from("grups-output");
        let mut fixtures: HashMap<&'static str, Fixture> = HashMap::new();

        // --- Parse module
        let Some(module) = self.module else {
            bail!("No module selected");
        };
        args.push(module.to_string());
        

        // --- pedigree-sims specific
        if module == "pedigree-sims" {
            let Some(pedigree) = self.pedigree else {
                bail!("No specified pedigree")
            };
            fixtures.insert("pedigree", Fixture::copy(pedigree));
            args.push(format!("--pedigree {}", fixtures.get("pedigree").unwrap()));

            let Some(data_dir) = self.data_dir else {
                bail!("No specified data directory")
            };
            fixtures.insert("data-dir", Fixture::copy(data_dir));
            args.push(format!("--data-dir {}", fixtures.get("data-dir").unwrap()));

            let Some(recomb_dir) = self.recomb_dir else {
                bail!("No specified recombination directory");
            };
            fixtures.insert("recomb-dir", Fixture::copy(recomb_dir));
            args.push(format!("--recomb-dir {}", fixtures.get("recomb-dir").unwrap()));

            if let Some(mode) = self.mode {
                let mode_str = match mode {
                    parser::Mode::Vcf     => "vcf",
                    parser::Mode::Fst     => "fst",
                    parser::Mode::FstMmap => "fst-mmap",
                };

                args.push(format!("--mode {mode_str}"));
            }

            match self.reps {
                Some(arg) => args.push(format!("--reps {arg}")),
                None      => args.push(format!("--reps {TEST_PEDIGREE_SIMS_REPS}"))
            }

            match self.seed {
                Some(arg) => args.push(format!("--seed {arg}")),
                None      => args.push(format!("--seed {TEST_PEDIGREE_SIMS_SEED}")),
            }

            if let Some(contam_pop) = self.contam_pop {
                args.push(format!("--contam-pop {contam_pop}"));
            }

            if self.sex_specific_mode {
                args.push("--sex-specific-mode".to_string());
            }
        }

        // ---- pwd-from-stdin or pedigree-sims
        if matches!(module, "pwd-from-stdin" | "pedigree-sims") {
            let Some(pileup) = self.pileup else {
                bail!("No specified pileup");
            };
            fixtures.insert("pileup", Fixture::copy(pileup));
            args.push(format!("--pileup {}", fixtures.get("pileup").unwrap()));
            filestem = Some(Path::new(pileup).file_stem().unwrap().to_string_lossy().to_string());

            if let Some(targets) = self.targets {
                fixtures.insert("targets", Fixture::copy(targets));
                args.push(format!("--targets {}", fixtures.get("targets").unwrap()));
            }

            if let Some(samples) = self.samples {
                args.push(format!("--samples {samples}"));
            }

            if self.exclude_transitions {
                args.push("--exclude-transitions".to_string());
            }

            if self.x_chromosome_mode {
                args.push("--x-chromosome-mode".to_string())
            }
        }

        if matches!(module, "pwd-from-stdin" | "pedigree-sims" | "fst") {
            if let Some(dir) = self.output_dir {
                output_dir = dir.to_string();
            }
            fixtures.insert("output-dir", Fixture::blank(&output_dir));
            args.push(format!("--output-dir {}", fixtures.get("output-dir").unwrap()));

            if self.overwrite {
                args.push("--overwrite".to_string());
            }
        }

        if matches!(module, "fst") {
            let Some(vcf_dir) = self.vcf_dir else {
                bail!("No specified vcf directory");
            };
            fixtures.insert("vcf-dir", Fixture::copy(vcf_dir));
            args.push(format!("--vcf-dir {}", fixtures.get("vcf-dir").unwrap()));

            if let Some(pop_subset) = self.pop_subset {
                args.push(format!("--pop-subset {pop_subset}"));
            }

            if self.compute_pop_afs {
                args.push("--compute-pop-afs".to_string());
            }
        }

        Ok(GrupsRunner {args, filestem, fixtures, module: module.to_string()})
    }

}


pub struct GrupsRunner {
    args: Vec<String>,
    filestem: Option<String>,
    fixtures: HashMap<&'static str, Fixture>,
    module: String,
}

impl GrupsRunner {
    pub fn run(&self) {
        let args = self.args.join(" ");
        let cli = parser::Cli::parse_from(args.split_whitespace());
        grups_rs::run(cli).expect("Failed to run grups using stringified CLI Args");
    }

    pub fn get_fixture(&self, key: &str) -> &Fixture {
        self.fixtures.get(key).unwrap()
    }

    fn filestem(&self) -> String {
        let output_dir = self.get_fixture("output-dir");
        let filestem = self.filestem.as_ref().unwrap();
        format!("{output_dir}/{filestem}")
    }

    #[allow(unused)]
    pub fn simulations_filestem(&self) -> String {
        let output_dir = self.get_fixture("output-dir");
        let filestem = self.filestem.as_ref().unwrap();
        format!("{output_dir}/simulations/{filestem}")

    }

    pub fn output_pwd(&self) -> Option<String> {
        let output = format!("{}.pwd", self.filestem());
        match self.module.as_str() {
            "pedigree-sims" => Some(output),
            _               => None
        }
    }

    pub fn output_results(&self) -> Option<String> {
        let output = format!("{}.result", self.filestem());
        match self.module.as_str() {
            "pwd-from-stdin" |"pedigree-sims" => Some(output),
            _                                 => None
        }
    }

    #[allow(unused)]
    pub fn output_probs(&self) -> Option<String> {
        let output = format!("{}.probs", self.filestem());
        match self.module.as_str() {
            "pwd-from-stdin" |"pedigree-sims" => Some(output),
            _                                 => None
        }
    }

    #[allow(unused)]
    pub fn args(&self) -> Vec<String> {
        self.args.clone()
    }
}
