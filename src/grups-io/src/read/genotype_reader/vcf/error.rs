use thiserror::Error;

#[derive(Error, Debug)]
pub enum VCFReaderError {
    #[error("Invalid or missing file extension. Accepted format are ['.vcf', '.vcf.gz']")]
    InvalidFileExt,

    #[error("Failed to open VCF file")]
    Open,

    #[error("Failed to parse VCF line.")]
    ParseLine, 

    #[error("Failed to read the contents of the VCF")]
    FillBuffer(#[source] std::io::Error),

    #[error("VCF Field contains invalid UTF8 data")]
    InvalidField(#[source] std::str::Utf8Error),

    #[error("Failed to jump to the next line")]
    SkipLineError(#[source] anyhow::Error),

    #[error("Failed to convert field into a valid chromosome index")]
    ParseChrError(#[source] genome::coordinate::ChrIdxError),

    #[error("Failed to convert field into a valid genomic position")]
    ParsePosError(#[source] genome::coordinate::ParsePositionError),

    #[error("Unexpected field index when performing operation. Must be {0}, was {1}")]
    InvalidFieldIdx(usize, usize),
}