use std::{fs::File, io::{Write, BufWriter}, path::Path};
use anyhow::Result;
use regex::Regex;

use located_error::LocatedError;

pub mod error;
pub use error::WriterError;

/// THE field separator used for this generic writer.
pub const WRITER_SEPARATOR: &str = "\t";

/// A generic file writer.
/// - source: Boxed `BufWriter` (can either handle file-writing, or stdout).
pub struct GenericWriter<'a> {
    source: BufWriter<Box<dyn Write + 'a>>
}

impl<'a> GenericWriter<'a>{
    /// Instantiate a new `Writer`, linked to a file.
    /// 
    /// # Errors
    /// if `path` is either an invalid file, or the user does not have the proper
    /// UNIX permissions to write at this location.
    pub fn new(path: Option<impl AsRef<Path>>) -> Result<GenericWriter<'a>>{
        use WriterError::{IOError};
        Ok(GenericWriter{ source: match path {
            Some(path) => {
                let file = File::create(path).map_err(IOError).loc("While creating file")?;
                BufWriter::new(Box::new(file))
            },
            None => {
                BufWriter::new(Box::new(std::io::stdout()))
            }
        }})
    }

    /// Write the contents of a generic iterator within a file/stdout.
    /// one Iteration step = one line.
    /// 
    /// # Behavior
    /// 
    /// For each item of the iterator, `write_iter` will search for the regular expression
    /// `[ ]+-[ ]+` and replace matches with `\t`. This effectively removes "Pretty-print" 
    /// from the output.
    /// 
    /// # Errors
    /// - If any of the Items within `iter` fails to get written within the file.
    /// 
    /// # Panics
    /// - if parsing the regex required to delete pretty-print characters fails.
    /// 
    pub fn write_iter<T, I>(&mut self, iter: T) -> Result<()>
    where   T: IntoIterator<Item = I>,
            I: std::fmt::Display,
    {
        // @TODO: Lazy static this regex!
        let re = Regex::new(r"[ ]+-[ ]+").unwrap(); // Remove pretty print trailing and leading whitespace.
        iter.into_iter()
            .map(|obj| self.source.write(re.replace_all(&format!("{}\n", obj), WRITER_SEPARATOR).as_bytes()))
            .collect::<Result<Vec<usize>, _>>()
            .map_err(WriterError::IOError)
            .loc("While writing contents into file")?;

        self.source.flush().loc("While flushing buffer contents of Writer")
    }
}
