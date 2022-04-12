use std::fs::File;
use std::io::{self, BufRead, Read};

use crate::genome::SNPCoord;

use std::collections::HashSet;
use std::error::Error;

#[derive(Debug)]
pub enum SNPReaderError {
    InvalidFileFormatError(String),
    FileNotFoundError(String, String),
}
impl Error for SNPReaderError {}

impl std::fmt::Display for SNPReaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SNPReaderError::InvalidFileFormatError(ext) => write!(f, "Cannot handle targetfile format: {}", ext),
            SNPReaderError::FileNotFoundError(path, err) => write!(f, "{}: {}",path, err)
        }
    }
}

// https://medium.com/bgpkit/write-generic-file-reader-in-rust-ad6408cb086a
// --> http::request::Request
pub struct SNPReader<'a> {
    source: Box<dyn BufRead + 'a>,
    columns: Vec<usize>,
    sep: String,
}

impl<'a> SNPReader<'a> {
    pub fn new(path: &str) -> Result<SNPReader<'a>, SNPReaderError> {
        let (columns, sep): (Vec<usize>, String) = Self::get_file_format(path)?;

        match File::open(path).map(|file| SNPReader {
            source: Box::new(io::BufReader::new(file)),
            columns,
            sep
        }) {
            Err(e) => Err(SNPReaderError::FileNotFoundError(path.to_string(), e.to_string())),
            Ok(reader) => Ok(reader)
        }
    }

    /// Extract file extension from the file path and return the appropriate columns and
    /// field-separator. 
    /// Supported formats: [.snp, .vcf, .txt, .csv, .tsv]
    /// 
    /// TODO: convert separators from String -> &str
    fn get_file_format(path: &str) -> Result<(Vec<usize>,String), SNPReaderError> {
        let file_type = path.split(".").collect::<Vec<&str>>().last().unwrap().clone();
        let output = match file_type {
            // ([CHR, POS, REF, ALT], SEP)
            "snp" => {(vec![1,3,4,5]," ".to_string())}
            "vcf" => {(vec![0,1,3,4],"\t".to_string())}
            "txt" => {(vec![0,1,2,3]," ".to_string())}
            "csv" => {(vec![0,1,2,3],",".to_string())}
            "tsv" => {(vec![0,1,2,3],"\t".to_string())}
            t => return Err(SNPReaderError::InvalidFileFormatError(t.to_string()))
        };
        Ok(output)
    }

    /// Read from a BufReader and convert each line to a HashMap of SNPCoord.
    /// Input: - path (string): path and filename to a snp coordinates file.
    ///
    /// Return: HashMap of structs 'SNPCoord'
    /// 
    /// TODO: This should take a &self reference.
    pub fn hash_target_positions(self) -> Result<HashSet<SNPCoord>, Box<dyn Error>> {
        let mut target_positions : HashSet<SNPCoord> = HashSet::new(); // Output
        for line in self.source.lines() {
            let line = line?;
            match &line[..1] {     // Skip comments.
                "#" => {continue}
                _   => {}
            }
            let split_line: Vec<&str>    = line.split(&self.sep).filter(|x| x!=&"").collect();
            let chromosome: u8           = split_line[self.columns[0]].parse()?;
            let position  : u32          = split_line[self.columns[1]].parse()?;
            let reference : Option<char> = split_line[self.columns[2]].parse().ok();
            let alternate : Option<char> = split_line[self.columns[3]].parse().ok();
            
            let coordinate: SNPCoord     = SNPCoord {chromosome, position, reference, alternate};
            target_positions.insert(coordinate);
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

pub struct Writer<'a> {
    source: BufWriter<Box<dyn Write + 'a>>
}

use std::io::Write;
use std::io::BufWriter;

impl<'a> Writer<'a>{
    pub fn new(path: Option<String>) -> Result<Writer<'a>, Box<dyn Error>>{
        Ok(Writer { source: match path {
            Some(path) => {
                //let file = File::create(path)?;
                BufWriter::new(Box::new(File::create(path)?))
            },
            None => {
                BufWriter::new(Box::new(std::io::stdout()))
            }
        }})
    }

    pub fn write_iter<T, I>(&mut self, iter: T) -> Result<(), std::io::Error>
    where
        T: IntoIterator<Item = I>,
        I: std::fmt::Display,
    {
        iter.into_iter()
            .map(|obj| self.source.write(format!("{}\n", obj).as_bytes()))
            .collect::<Result<Vec<usize>, std::io::Error>>()?;
        self.source.flush()
    }
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