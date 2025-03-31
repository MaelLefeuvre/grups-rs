use std::{fs::File, io::{Write, BufWriter}, path::Path};
use anyhow::Result;
use regex::Regex;
use lazy_static::lazy_static;

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
        use WriterError::IOError;
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
        // Remove pretty print trailing and leading whitespace
        lazy_static! {
            static ref RE: Regex = Regex::new(r"[ ]+-[ ]+").expect("Failed to parse regex.");
        }
        iter.into_iter()
            .map(|obj| self.source.write(RE.replace_all(&format!("{obj}\n"), WRITER_SEPARATOR).as_bytes()))
            .collect::<Result<Vec<usize>, _>>()
            .map_err(WriterError::IOError)
            .loc("While writing contents into file")?;

        self.source.flush().loc("While flushing buffer contents of Writer")
    }
}

#[cfg(test)]
mod tests {
    use genome::coordinate::Coordinate;

    use super::*;
    #[test]
    fn write_file() -> anyhow::Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.txt");
        let mut writer = GenericWriter::new(Some(&path))?;

        let test_vec = vec![Coordinate::new(10, 10000)];
        writer.write_iter(&test_vec)?;

        let got = std::io::read_to_string(File::open(path)?)?;
        assert_eq!(got.replace('\n', ""), test_vec[0].to_string());
        Ok(())
    }
}