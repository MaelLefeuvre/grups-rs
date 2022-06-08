#[derive(Debug, Clone)]
pub struct PedigreeParams {
    pub snp_downsampling_rate : f64,
    pub af_downsampling_rate  : f64,
    pub seq_error_rate        : Option<[f64; 2]>,
    pub contam_rate           : [f64; 2],
}

impl PedigreeParams {
    pub fn new(snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: Option<[f64; 2]>, contam_rate: [f64; 2]) -> Self {
        PedigreeParams{snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate}
    }
}