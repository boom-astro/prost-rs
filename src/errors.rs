use thiserror::Error;

#[derive(Error, Debug)]
pub enum ProstError {
    #[error("no candidates found within search radius")]
    NoCandidates,

    #[error("invalid galaxy shape parameters: {0}")]
    InvalidShape(String),

    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    #[error("numerical error: {0}")]
    NumericalError(String),

    #[error("cutout retrieval failed: {0}")]
    CutoutError(String),

    #[error("source extraction failed: {0}")]
    SourceExtractionError(String),

    #[error("morphology fitting failed: {0}")]
    MorphologyError(String),
}
