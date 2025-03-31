use std::{fs::File, io::{self, BufRead, BufReader, Read}, str::FromStr};
use genome::{coordinate::ChrIdx, snp::Allele, SNPCoord};
use located_error::*;
use anyhow::Result;
use log::info;

use ahash::AHashSet;

pub mod error;
pub use error::SNPReaderError;

pub const SNPREADER_VALID_FILE_FORMATS: [&str; 5] = ["snp", "vcf", "txt", "csv", "tsv"];

#[derive(PartialEq, Eq, Debug)]
pub enum SNPReaderMode {Snp, Vcf, Txt, Tsv, Csv}

impl FromStr for SNPReaderMode {
    type Err = SNPReaderError;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "snp" => Ok(Self::Snp),
            "vcf" => Ok(Self::Vcf),
            "txt" => Ok(Self::Txt),
            "csv" => Ok(Self::Csv),
            "tsv" => Ok(Self::Tsv),
            other => Err(Self::Err::InvalidFileFormat(other.to_string()))
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
    mode: SNPReaderMode,
    columns: [usize; 4],
    sep: char
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
        let (mode, columns, sep) = Self::get_file_format(path)?;

        File::open(path)
            .map(|file| SNPReader {source: Box::new(BufReader::new(file)), mode, columns, sep})
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
    fn get_file_format(path: &str) -> Result<(SNPReaderMode, [usize; 4], char)> {
        let file_type: &str = path.split('.')
            .collect::<Vec<&str>>()
            .last()
            .ok_or(SNPReaderError::MissingExtension)
            .loc("While attempting to extract file extension")?;
        let mode = SNPReaderMode::from_str(file_type).loc("While attempting to retrieve input file format")?;
        let (columns, sep) = match mode { // ([CHR, POS, REF, ALT], SEP)
            SNPReaderMode::Snp => ([1,3,4,5], ' ' ),
            SNPReaderMode::Vcf => ([0,1,3,4], '\t'),
            SNPReaderMode::Txt => ([0,1,2,3], ' ' ),
            SNPReaderMode::Csv => ([0,1,2,3], ',' ),
            SNPReaderMode::Tsv => ([0,1,2,3], '\t'),
        };
        Ok((mode, columns, sep))
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
            let mut chromosome: ChrIdx = split_line[self.columns[0]].parse().with_loc(context)?;
            let position  : u32        = split_line[self.columns[1]].parse().with_loc(context)?;
            let reference : Allele     = split_line[self.columns[2]].parse().with_loc(context)?;
            let alternate : Allele     = split_line[self.columns[3]].parse().with_loc(context)?;

            if chromosome == ChrIdx(23) && self.mode == SNPReaderMode::Snp {
                chromosome = ChrIdx(b'X')
            };
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
    use std::io::Write;

    macro_rules! check_file_format {
        ($cols:expr, $sep:expr, $ext:expr) => {
            let path = format!("./targets /v50.0_1240K_public.{}", $ext);
            assert_eq!((SNPReaderMode::from_str($ext).unwrap(), $cols, $sep), SNPReader::get_file_format(&path).expect("Failed to obtain file format"));
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

    #[test]
    fn test_hash_position_tsv() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.tsv");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
            14\t1565489\tA\tC\n\
            15\t1500000\tG\tT\n\
            chrX\t17000\tG\tA\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(false)?;

        assert!(positions.contains(&SNPCoord::new(14, 1565489, Allele::A, Allele::C)));
        assert!(positions.contains(&SNPCoord::new(15, 1500000, Allele::G, Allele::T)));
        assert!(positions.contains(&SNPCoord::new(b'X', 17000, Allele::G, Allele::A)));
        Ok(())
    }

    #[test]
    fn test_hash_position_txt() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.txt");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
            14 1565489 A C\n\
            15 1500000 G T\n\
            chrX 15000 A G\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(false)?;

        assert!(positions.contains(&SNPCoord::new(14, 1565489, Allele::A, Allele::C)));
        assert!(positions.contains(&SNPCoord::new(15, 1500000, Allele::G, Allele::T)));
        assert!(positions.contains(&SNPCoord::new(b'X', 15000, Allele::A, Allele::G)));
        Ok(())
    }

    #[test]
    fn test_hash_position_csv() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.csv");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
            14,1565489,A,C\n\
            15,1500000,G,T\n\
            chrX,18000,A,C\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(false)?;

        assert!(positions.contains(&SNPCoord::new(14, 1565489, Allele::A, Allele::C)));
        assert!(positions.contains(&SNPCoord::new(15, 1500000, Allele::G, Allele::T)));
        assert!(positions.contains(&SNPCoord::new(b'X', 18000, Allele::A, Allele::C)));
        Ok(())
    }

    #[test]
    fn test_hash_position_vcf() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.vcf");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
            ##fileformat=VCFv4.1\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\n\
            15\t60026\t.\tA\tC\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.5;AFR_AF=0.5;EUR_AF=0.5;SAS_AF=0.5;EAS_AF=0.5;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\n\
            16\t60057\t.\tG\tT\t100\tPASS\tAC=1582;AF=0.31;AMR_AF=0.02;AFR_AF=0.04;EUR_AF=0.06;SAS_AF=0.08;EAS_AF=0.010;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\n\
            X\t60083\t.\tG\tA\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.43;AFR_AF=0.38;EUR_AF=0.99;SAS_AF=0.61;EAS_AF=0.01;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(false)?;

        assert!(positions.contains(&SNPCoord::new(15,   60026, Allele::A, Allele::C)));
        assert!(positions.contains(&SNPCoord::new(16,   60057, Allele::G, Allele::T)));
        assert!(positions.contains(&SNPCoord::new(b'X', 60083, Allele::G, Allele::A)));
        Ok(())
    }

    #[test]
    fn test_hash_position_snp() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.snp");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
           rs3094315     1        0.020130          752566 G A\n\
         rs115861570    16        0.000623          155344 G A\n\
      snp_23_2779345    23        0.211432         2779345 C T\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(false)?;

        assert!(positions.contains(&SNPCoord::new(1, 752566, Allele::G, Allele::A)));
        assert!(positions.contains(&SNPCoord::new(16, 155344, Allele::G, Allele::A)));
        assert!(positions.contains(&SNPCoord::new(88, 2779345, Allele::C, Allele::T)));
        Ok(())
    }

    #[test]
    fn test_hash_position_transitions() -> Result<()> {
        let tmpdir        = tempfile::tempdir()?;
        let path          = tmpdir.path().join("targets.csv");
        let provided_path = path.to_str().expect("Invalid path");
        let mut file      = File::create(&path)?;
        writeln!(file, "\
            1,1,A,C\n\
            1,2,A,G\n\
            1,3,A,T\n\
            1,4,C,A\n\
            1,5,C,G\n\
            1,6,C,T\n\
            1,7,G,A\n\
            1,8,G,C\n\
            1,9,G,T\n\
            1,10,T,A\n\
            1,11,T,C\n\
            1,12,T,G\
        ")?;
        let reader        = SNPReader::new(provided_path)?;
        let positions = reader.hash_target_positions(true)?;
        let transitions_pos =  [2, 6, 7, 11];
        for position in 1..=12 {
            let snpcoord = SNPCoord::new(1, position, Allele::N, Allele::N);
            let coord_was_retained = positions.contains(&snpcoord); 
            match transitions_pos.contains(&position) {
                true  => assert!(!coord_was_retained),
                false => assert!(coord_was_retained),
            }
        }
        Ok(())
    }
}