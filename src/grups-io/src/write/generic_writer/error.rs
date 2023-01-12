use thiserror::Error;

#[derive(Error, Debug)]
pub enum WriterError {
    #[error("Failed to write to file: inner writter return io error")]
    IOError(#[from] std::io::Error)
}