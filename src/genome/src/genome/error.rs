
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FastaIndexReaderError {
    #[error("{path}: {err}")]
    FileNotFound{path: String, err: String},

    #[error("{path}: Line {line} - Invalid line format [{err}]")]
    ParseInt{path: String, line: u8, err: String},

    #[error("Invalid fasta file format:\n - expected [.fa |.fa.fai |.fasta | .fasta.fai]\n - got '.{ext}'")]
    InvalidExt{ext: String},

    #[error("At line {idx}: Failed to parse fasta index line - got [{err}]")]
    ParseLine{idx: usize, err: std::io::Error}
}

