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

}