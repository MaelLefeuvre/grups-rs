#![allow(unused)]
use pwd_from_stdin::comparisons::Comparisons as PileupComparisons;

#[derive(Default)]
struct PedigreeSimulatorBuilder<'a> {
    // A Reader: This reader must be able to change from one file to another
    //            The reader keeps track of the list of files it should iterate on.
    
    // A HashMap of pedigrees

    // A Genetic Map
    // A HashMap of Previous_positions

    // an RNG ?

    // Files
    data_directory       : Option<String>,
    recombination_map    : Option<String>,
    pileup_comparisons   : Option<&'a PileupComparisons>,
    panel_definition_file: Option<String>,

    // Pedigree
    pedigree             : Option<String>,
    reps                 : Option<usize>,
    // Population definition
    pedigree_pop         : Option<String>,
    contam_pops          : Option<String>,
    contam_num_inds      : Option<Vec<usize>>,

    // Params
    snp_downsampling_rate: Option<f64>,
    af_downsampling_rate : Option<f64>,
    seq_error_rates      : Option<Vec<Vec<f64>>>,
    contam_rates         : Option<Vec<Vec<f64>>>,    
}

impl PedigreeSimulatorBuilder<'_> {
    fn new() -> Self {
        Self::default()
    }

    fn pedigree_pop(&mut self, pop: &str) -> &mut Self {
        self.pedigree_pop = Some(pop.to_string());
        self
    }

    fn pedigree_reps(&mut self, n: usize) -> &mut Self {
        self.reps = Some(n);
        self
    }

    fn contam_pops(&mut self, contam_pops: &str) -> &mut Self {
        self.contam_pops = Some(contam_pops.to_string());
        self
    }

    fn contam_num_inds(&mut self, n: Vec<usize>) -> &mut Self {
        self.contam_num_inds = Some(n);
        self
    }
    fn snp_downsampling_rate(&mut self, rate: f64 ) -> &mut Self {
        self.snp_downsampling_rate = Some(rate);
        self
    }

    fn af_downsampling_rate(&mut self, rate: f64) -> &mut Self {
        self.af_downsampling_rate = Some(rate);
        self
    }

    fn seq_error_rates(&mut self, rate: Vec<Vec<f64>>) -> &mut Self {
        self.seq_error_rates = Some(rate);
        self
    }

    fn contam_rates(&mut self, rates: Vec<Vec<f64>>) -> &mut Self {
        self.contam_rates = Some(rates);
        self
    }
}