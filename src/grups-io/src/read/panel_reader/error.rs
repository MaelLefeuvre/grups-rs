use std::path::PathBuf;

use thiserror::Error;

const COMMON_MSG: &str = "Please place a valid '.panel' file within the specified '--data-dir', or directly specify the file using '--panel'";

#[derive(Error, Debug)]
pub enum PanelReaderError {
    #[error("Could not find a candidate panel definition file within the requested directory. {}", COMMON_MSG)]
    NotFound,

    #[error("Found multiple candidate Panel definition files: {0:#?}. Please specify the relevant file using '--panel'")]
    MultipleFound(Vec<PathBuf>),

    #[error("Failed to retrieve a valid random SampleTag from our panel, using population tag {0}. Note that this error can happen when providing '--contam-pop' with a population identifier that is invalid and/or missing from the panel definition file.")]
    MissingContaminant(String),

    #[error("Could not fetch random sample using population tag: '{0}'")]
    MissingSample(String),

    #[error("Exhausted the number of available samples to populate the pedigree")]
    ExhaustedPanel,

    #[error("Line {0} of the provided .panel definition file does not contain the appropriate number of fields. Each line of this file should contain 3 to 4 fields: <ID (required)> <super-population-tag (required)> <population tag (required)> <sex (optional)> (e.g. 'HG00096 GBR EUR male')")]
    InvalidNumberOfFields(usize),

    #[error("Line {0} - Field 1 of the provided .panel definition file carries an invalid Individual-ID.")]
    ParseIndividualId(usize),

    #[error("Line {0} - Field 2 of the provided .panel definition file carries an invalid Superpopulation ID.")]
    ParseSuperPop(usize),

    #[error("Line {0} - Field 3 of the provided .panel definition file carries an invalid Population ID.")]
    ParsePop(usize),

}