use std::path::Path;
use std::fs::File;
use std::io::{self, BufRead, Read};

use genome::SNPCoord;

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io::Write;
use std::io::BufWriter;
use std::path::PathBuf;

use log::{info, trace};

#[derive(Debug)]
pub enum SNPReaderError {
    InvalidFileFormat(String, Vec<&'static str>),
    FileNotFound(String, String),
    MissingExtension(Vec<&'static str>),
    MissingAltRef(SNPCoord),
}
impl Error for SNPReaderError {}

impl std::fmt::Display for SNPReaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::InvalidFileFormat(ext, file_formats) => {
                write!(f, "Cannot handle target file format: {ext} \
                    Please provide '--targets' with either one of the accepted file formats: {file_formats:?}"
                )
            },
            Self::FileNotFound(path, err) => {
                write!(f, "{}: {}", path, err)
            },
            Self::MissingExtension(file_formats) => {
                write!(f, "The provided targets file is missing a file extension. \
                    Please provide '--targets' with a file ending with one of the accepted file formats: \n
                    {file_formats:?}"
                )
            },
            Self::MissingAltRef(nuc) => {
                write!(f, "Target SNP Coordinate [{} {}] is missing its reference and/or alternate allele \
                    information when it is required.",
                    nuc.chromosome, nuc.position
                )
            }


        }
    }
}


// https://medium.com/bgpkit/write-generic-file-reader-in-rust-ad6408cb086a
// --> http::request::Request
/// Generic file reader for snp coordinates files. 
/// - source  : `BufReader`. Currently Boxed, because it might be interesting to incorporate HTTP requests.
/// - columns : a vector of size 4. indidicating the column indices that we're targeting (chr, pos, ref, alt)
/// - sep     : the expected field separator for the corresponding file format.
/// 
/// ## Accepted file formats: 
///    EXT  NAME                 CHR  POS  REF  ALT  SEP
/// - .snp  eigenstrat format    1    3    4    5    " "
/// - .vcf  variant call format  0    1    3    4    "\t"
/// - .csv  comma separated      0    1    2    3    "," 
/// - .tsv  tab separated        0    1    2    3    "\t"
/// - .txt  space separated      0    1    2    3    " "
/// 
/// # Traits : `Read` and `BufRead`.
pub struct SNPReader<'a> {
    source: Box<dyn BufRead + 'a>,
    columns: Vec<usize>,
    sep: String,
}

impl<'a> SNPReader<'a> {
    /// Instantiate a new `SNPReader`
    /// 
    /// # Errors
    /// - if `path` targets an invalid location (`FileNotFound` or `PermissionDenied`)
    pub fn new(path: &str) -> Result<SNPReader<'a>, SNPReaderError> {
        use SNPReaderError::FileNotFound;
        let (columns, sep): (Vec<usize>, String) = Self::get_file_format(path)?;

        match File::open(path).map(|file| SNPReader {
            source: Box::new(io::BufReader::new(file)),
            columns,
            sep
        }) {
            Err(e) => Err(FileNotFound(path.to_string(), e.to_string())),
            Ok(reader) => Ok(reader)
        }
    }

    /// Extract file extension from the file path and return the appropriate columns and
    /// field-separator. 
    /// Supported formats: [.snp, .vcf, .txt, .csv, .tsv]
    /// 
    /// # Errors
    ///  - `InvalidFileFormat` if `path` is not carrying a supported file extension.
    ///  - `MissingExtension` if `path` does not contain a file extension
    /// @TODO: convert separators from String -> &str
    fn get_file_format(path: &str) -> Result<(Vec<usize>,String), SNPReaderError> {
        use SNPReaderError::{InvalidFileFormat, MissingExtension};
        let accepted_file_formats = ["snp", "vcf", "txt", "csv", "tsv"];
        let file_type: &str = path.split('.')
            .collect::<Vec<&str>>()
            .last()
            .ok_or_else(|| MissingExtension(accepted_file_formats.to_vec()))?;
        let output = match file_type {
            // ([CHR, POS, REF, ALT], SEP)
            "snp" => (vec![1,3,4,5], " ".to_string() ),
            "vcf" => (vec![0,1,3,4], "\t".to_string()),
            "txt" => (vec![0,1,2,3], " ".to_string() ),
            "csv" => (vec![0,1,2,3], ",".to_string() ),
            "tsv" => (vec![0,1,2,3], "\t".to_string()),
            t => return Err(InvalidFileFormat(t.to_string(), accepted_file_formats.to_vec()))
        };
        Ok(output)
    }

    /// Read from a  `BufReader` and convert each line to a `HashMap` of `SNPCoord`.
    /// Input: - path (string): path and filename to a snp coordinates file.
    /// Return: `HashMap` of structs `SNPCoord`
    /// 
    /// # Errors
    /// - Can return `ParseIntError` or `ParseCharError` if any of the fields contains invalid data.
    /// - Returns an error if any line fails to get parsed into a string.
    /// 
    /// # Panics
    /// - If `exclude_transitions` is on, and `self.source` does not contain any REF/ALT information.
    /// 
    /// # @TODO: 
    ///  - This should take a &self reference.
    pub fn hash_target_positions(self, exclude_transitions: bool) -> Result<HashSet<SNPCoord>, Box<dyn Error>> {
        let mut target_positions : HashSet<SNPCoord> = HashSet::new(); // Output
        for line in self.source.lines() {
            let line = line?;
            if &line[..1] == "#" {continue} // Skip comment
            let split_line: Vec<&str>    = line.split(&self.sep).filter(|x| x!=&"").collect();
            let chromosome: u8           = split_line[self.columns[0]].parse()?;
            let position  : u32          = split_line[self.columns[1]].parse()?;
            let reference : Option<char> = split_line[self.columns[2]].parse().ok();
            let alternate : Option<char> = split_line[self.columns[3]].parse().ok();
            
            let coordinate: SNPCoord     = SNPCoord {chromosome, position, reference, alternate};

            if exclude_transitions && !coordinate.has_known_alleles() {
                return Err("Cannot filter-out transitions if reference and alternate alleles is unknown. \
                Please provide '--targets' with a file containing known reference and alternate alleles.".into())
            }

            target_positions.insert(coordinate);
        }

        let transitions = vec![['A', 'G'], ['G', 'A'], ['C', 'T'], ['T', 'C']];
        if exclude_transitions {
            info!("Filtering transitions from targets file.");
            let before = target_positions.len();
            target_positions.retain(|nuc| {
                !transitions.contains(&[
                    nuc.reference.unwrap(), 
                    nuc.alternate.unwrap()
                ])
            });
            let after = target_positions.len();
            info!("{} transitions filtered out. (Before: {} - After: {})", before-after, before, after);
        }

        Ok(target_positions)
    }

}

impl<'a> Read for SNPReader<'a> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.source.read(buf)
    }
}

impl<'a> BufRead for SNPReader<'a> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.source.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.source.consume(amt);
    }
}

/// A generic file writer.
/// - source: Boxed `BufWriter` (can either handle file-writing, or stdout).
pub struct Writer<'a> {
    source: BufWriter<Box<dyn Write + 'a>>
}

impl<'a> Writer<'a>{
    /// Instantiate a new `Writer`, linked to a file.
    /// 
    /// # Errors
    /// if `path` is either an invalid file, or the user does not have the proper
    /// UNIX permissions to write at this location.
    pub fn new(path: Option<String>) -> Result<Writer<'a>, Box<dyn Error>>{
        Ok(Writer { source: match path {
            Some(path) => {
                BufWriter::new(Box::new(File::create(path)?))
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
    pub fn write_iter<T, I>(&mut self, iter: T) -> Result<(), std::io::Error>
    where
        T: IntoIterator<Item = I>,
        I: std::fmt::Display,
    {
        const WRITER_SEPARATOR: &str = "\t";
        let re = regex::Regex::new(r"[ ]+-[ ]+").unwrap(); // Remove pretty print trailing and leading whitespace.
        iter.into_iter()
            .map(|obj| self.source.write(re.replace_all(format!("{}\n", obj).as_str(), WRITER_SEPARATOR).as_bytes()))
            .collect::<Result<Vec<usize>, std::io::Error>>()?;
        self.source.flush()
    }
}


/// Obtain predefined filenames from the given output directory and pileup file. Return a `HashMap` with 
/// K: file-ext, V: Filepath. Right now this only outputs a single file, but is easily scalable.
///  - {out-dir}/{file-prefix}.pwd -> where summary statistics are printed for pairwise differences. 
/// 
/// # Errors
///  - If creating the parent directories of `file_prefix` is required: will throw a `PermissionDenied`
///    if the user does not have the proper UNIX permissions
/// 
/// # Panics
/// - If failing to convert `file_prefix` from `PathBuf` to `&str`
/// 
pub fn get_results_file_prefix<'a>(file_prefix: &'a mut PathBuf, file_ext: Vec<&'a str>) -> std::io::Result<HashMap<&'a str, String>> {
    // Create output directory. Early return if we can't create it
    std::fs::create_dir_all(file_prefix.parent().unwrap_or(file_prefix))?; 

    // Generate a HashMap of filepaths from the file_prefix
    let mut outfiles_hash = HashMap::new();
    for ext in file_ext {
        file_prefix.set_extension(ext);
        trace!("Output File: {file_prefix:?}");
        outfiles_hash.insert(ext, String::from(file_prefix.to_str().unwrap()));
    }
    Ok(outfiles_hash)
}

/// Simple enum for `get_output_files()`
///  Suffix => `HashMap` will use provided file suffixes as keys
///  Key    => `HashMap` will use provided file extensions as keys 
#[derive(Clone, Copy)]
pub enum FileKey{Suffix, Ext}

/// Obtain predefined filenames for jackknife blocks from the given output directory and pileup file. 
/// Return a `HashMap` with (K (pair): "{ind1}-{ind2}", V (File): "{out dir}/blocks/{file prefix}.block"
/// ==> A new file is generated for each pair of individuals. 
/// 
/// # Errors
/// 
/// If creating the parent directory of `file_prefix` is required, will throw a `PermissionDenied` if
/// the user does not have the proper UNIX permissions 
/// 
/// # Panics
/// 
/// when failing to convert `file_prefix` from `PathBuf` to `&str`
/// 
pub fn get_output_files<'a>(file_prefix: &'a mut PathBuf, allow_overwrite: bool, sort: FileKey, suffixes: &[String], file_ext: &[&'a str]) -> std::io::Result<HashMap<String, String>> {
    // Create output directory. Early return if we can't create it
    let root_dir = file_prefix.parent().unwrap_or(file_prefix);
    match std::fs::create_dir_all(root_dir){
        Ok(()) => (),
        Err(err) => return Err(
            std::io::Error::new(std::io::ErrorKind::PermissionDenied,
            format!("Failed to create parent directory: {}. Got: {}", root_dir.to_str().unwrap_or_default(), err)
            )
        )
    }; 
    // Generate a HashMap of filepaths from the file_prefix
    let mut outfiles_hash = HashMap::new();
    for suffix in suffixes {
        let mut file = PathBuf::from(
            format!("{}{}{}.", // final dot is to fake an extension
                file_prefix.to_str().unwrap(),
                if suffix.is_empty() {""} else {"-"},
                suffix
            )
        );

        trace!("File: {file:?}");
        for ext in file_ext {
            file.set_extension(ext);
            can_write_file(allow_overwrite, &file)?;
            let key = match &sort {
                FileKey::Suffix => suffix.clone(),
                FileKey::Ext    => (*ext).to_string()
            };

            trace!("Output file: {file:?}");
            outfiles_hash.insert(key, String::from(file.to_str().unwrap()));
        }
    }
    Ok(outfiles_hash)
}

/// Check if a given file already exists ; raise an error if such is the case, and the user did not explicitly 
/// allow file overwriting.
/// # Errors
/// - If the provided `pathbuf` already exists and the user did not specifically allow for file
///   overwrite using the `--overwrite` argument
/// 
/// # Panics
/// - if the provided `pathbuf` fails to get parsed as a string.
/// 
/// # @TODO:
/// - This method is a duplicate of `parser::can_write_file`. I/O functions, methods and structs should be contained
///   is their own package.
pub fn can_write_file(overwrite: bool, pathbuf: &Path) -> std::io::Result<bool> {
    if ! overwrite && pathbuf.exists() {   // Check if this file already exists and/or if overwrite is allowed.
        return Err(std::io::Error::new(std::io::ErrorKind::AlreadyExists,
            format!("{:?} exists. use --overwrite to force.",
            pathbuf.to_str().unwrap()))
        )
    }
    Ok(true)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_reader_file_ext_snp(){
        let path = "./targets /v50.0_1240K_public.snp";
        assert_eq!((vec![1,3,4,5]," ".to_string()), SNPReader::get_file_format(&path).unwrap())
    }
    #[test]
    fn snp_reader_file_ext_vcf(){
        let path = "./targets /v50.0_1240K_public.vcf";
        assert_eq!((vec![0,1,3,4],"\t".to_string()), SNPReader::get_file_format(&path).unwrap())
    }
    #[test]
    fn snp_reader_file_ext_txt(){
        let path = "./targets /v50.0_1240K_public.txt";
        assert_eq!((vec![0,1,2,3]," ".to_string()), SNPReader::get_file_format(&path).unwrap())
    }

    #[test]
    fn snp_reader_file_ext_tsv(){
        let path = "./targets /v50.0_1240K_public.tsv";
        assert_eq!((vec![0,1,2,3],"\t".to_string()), SNPReader::get_file_format(&path).unwrap())
    }

    #[test]
    fn snp_reader_file_ext_csv(){
        let path = "./targets /v50.0_1240K_public.csv";
        assert_eq!((vec![0,1,2,3],",".to_string()), SNPReader::get_file_format(&path).unwrap())
    }
}