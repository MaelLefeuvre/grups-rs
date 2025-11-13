/// Wrapper struct containing constant pedigree parameters used during pedigree simulations.
/// # Fields:
/// - `snp_downsampling_rate`: probability of ignoring an SNP position during simulations.
/// 
/// - `af_downsampling_rate` : probability of fixating an SNP position during simulations. 
/// 
/// - `seq_error_rate`       : size-two array of probability of simulating a sequencing error during simulations. (seq_error_rate\[i\] corresponds to Individual\[i\]).
///   This field is optional, depending on whether or not the user provided a value.
///   - when `None`: The user did not provide any set value.seq_error_rate is computed at each position, using the corresponding Pileup Phred-score.
///   - when `Some`: The user has requested for a set sequencing error rate, across all SNP positions
/// 
/// - `contam_rate`          : probability of simulating a modern human contamination during simulations (contam_rate\[i\] corresponds to Individual\[i\])
/// 
#[allow(clippy::struct_field_names)]
#[derive(Debug, Clone)]
pub struct PedigreeParams {
    pub snp_downsampling_rate : f64,
    pub af_downsampling_rate  : f64,
    pub seq_error_rate        : Option<[f64; 2]>,
    pub contam_rate           : [f64; 2],
}

impl PedigreeParams {
    /// Instantiate a new PedigreeParams wrapper struct, from the user-provided parameters.
    pub fn new(snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: Option<[f64; 2]>, contam_rate: [f64; 2]) -> Self {
        PedigreeParams{snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate}
    }
}