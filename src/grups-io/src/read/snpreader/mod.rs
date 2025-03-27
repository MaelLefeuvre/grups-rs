use std::{fs::File, io::{self, Read, BufRead, BufReader}};
use genome::{coordinate::ChrIdx, snp::Allele, SNPCoord};
use located_error::*;
use anyhow::Result;
use log::info;

use ahash::AHashSet;

pub mod error;
pub use error::SNPReaderError;

pub const SNPREADER_VALID_FILE_FORMATS: [&str; 5] = ["snp", "vcf", "txt", "csv", "tsv"];


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
    columns: [usize; 4],
    sep: char,
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


impl<'a> SNPReader<'a> {
    /// Instantiate a new `SNPReader`
    /// 
    /// # Errors
    /// - if `path` targets an invalid location (`FileNotFound` or `PermissionDenied`)
    pub fn new(path: &str) -> Result<SNPReader<'a>> {
        use SNPReaderError::OpenFile;
        let (columns, sep): ([usize; 4], char) = Self::get_file_format(path)?;

        File::open(path)
            .map(|file| SNPReader {source: Box::new(BufReader::new(file)), columns, sep})
            .map_err(|e| OpenFile(path.to_string(), e))
            .loc("While attempting to create a new SNPReader")
    }

    /// Extract file extension from the file path and return the appropriate columns and
    /// field-separator. 
    /// Supported formats: [.snp, .vcf, .txt, .csv, .tsv]
    /// 
    /// # Errors
    ///  - `InvalidFileFormat` if `path` is not carrying a supported file extension.
    ///  - `MissingExtension` if `path` does not contain a file extension
    /// @TODO: convert separators from String -> &str
    fn get_file_format(path: &str) -> Result<([usize; 4], char)> {
        use SNPReaderError::{InvalidFileFormat, MissingExtension};
        let file_type: &str = path.split('.')
            .collect::<Vec<&str>>()
            .last()
            .ok_or(MissingExtension)
            .loc("While attempting to extract file extension")?;
        let output = match file_type {
                     // ([CHR, POS, REF, ALT], SEP)
            "snp" => ([1,3,4,5], ' ' ),
            "vcf" => ([0,1,3,4], '\t'),
            "txt" => ([0,1,2,3], ' ' ),
            "csv" => ([0,1,2,3], ',' ),
            "tsv" => ([0,1,2,3], '\t'),
            t => return Err(InvalidFileFormat(t.to_string()))
                .loc("While attempting to retrieve input file format")
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
    pub fn hash_target_positions(self, exclude_transitions: bool) -> Result<AHashSet<SNPCoord>> {
        let mut target_positions : AHashSet<SNPCoord> = AHashSet::new(); // Output
        let context = || "While hashing target positions";
        for line in self.source.lines() {
            let line = line?;
            if &line[..1] == "#" {continue} // Skip comment
            let split_line: Vec<&str> = line.split(self.sep).filter(|x| x!=&"").collect();
            let chromosome: ChrIdx    = split_line[self.columns[0]].parse().with_loc(context)?;
            let position  : u32       = split_line[self.columns[1]].parse().with_loc(context)?;
            let reference : Allele    = split_line[self.columns[2]].parse().with_loc(context)?;
            let alternate : Allele    = split_line[self.columns[3]].parse().with_loc(context)?;
            
            let coordinate: SNPCoord     = SNPCoord::try_new(chromosome, position, reference, alternate)?;

            if exclude_transitions && !coordinate.has_known_alleles() {
                return Err(SNPReaderError::UnknownAlleles).with_loc(context)
            }

            target_positions.insert(coordinate);
        }

        use genome::snp::Allele::*;
        let transitions = [[A, G], [G, A], [C, T], [T, C]];
        if exclude_transitions {
            info!("Filtering transitions from targets file.");
            let before = target_positions.len();
            target_positions.retain(|nuc| {
                !transitions.contains(&[nuc.reference, nuc.alternate])
            });
            let after = target_positions.len();
            info!("{} transitions filtered out. (Before: {} - After: {})", before - after, before, after);
        }

        Ok(target_positions)
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! check_file_format {
        ($cols:expr, $sep:expr, $ext:expr) => {
            let path = format!("./targets /v50.0_1240K_public.{}", $ext);
            assert_eq!(($cols, $sep), SNPReader::get_file_format(&path).expect("Failed to obtain file format"));
        }
    }

    #[test]
    fn snp_reader_file_ext_snp(){
        check_file_format!([1,3,4,5], ' ', "snp");
    }
    #[test]
    fn snp_reader_file_ext_vcf(){
        check_file_format!([0,1,3,4], '\t', "vcf");

    }
    #[test]
    fn snp_reader_file_ext_txt(){
        check_file_format!([0,1,2,3], ' ', "txt");
    }

    #[test]
    fn snp_reader_file_ext_tsv(){
        check_file_format!([0,1,2,3], '\t', "tsv");
    }

    #[test]
    fn snp_reader_file_ext_csv(){
        check_file_format!([0,1,2,3], ',', "csv");
    }
}