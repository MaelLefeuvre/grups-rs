use thiserror::Error;


#[derive(Error, Debug)]
pub enum SVMError {
    #[error("Failed to predict probabilities from the provided value(s).")]
    PredictProbs(#[source] libsvm::Error),

    #[error("Failed to predict boolean class from the provided value(s)")]
    PredictBool(#[source] libsvm::Error),

    #[error("Invalid SvmPredictor labels.")]
    InvalidLabels,
}

#[derive(Error, Debug)]
pub enum SvmBuilderError {

    #[error("Missing labels while attempting to build SVM")]
    MissingLabels,

    #[error("Missing features while attempting to build SVM")]
    MissingFeatures,

    #[error("Failed to fit SVM")]
    FitSvm(#[source] libsvm::Error),
}