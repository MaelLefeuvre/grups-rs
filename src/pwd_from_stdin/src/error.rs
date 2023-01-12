use thiserror::Error;

#[derive(Error, Debug)]
pub enum PwdFromStdinError {
    #[error("The use of '--exclude-transitions' requires an input SNP targets file, with known <REF> and <ALT> columns. Please provide this file, using the '--targets' argument")]
    MissingTargetPositions,

    #[error("Cannot filter known variants when REF/ALT allele are unknown! Please use a different file format.")]
    MissingKnownVariant
}