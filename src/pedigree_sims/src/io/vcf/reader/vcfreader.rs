use std::{
    io::{BufRead, BufReader, Read},
    path::Path,
    fs::File,
    error::Error,
     ops::Deref,
};


use crate::{
    io::{
        genotype_reader::GenotypeReader,
        vcf::sampletag::SampleTag,
    },
};
use gzp::{deflate::Bgzf, par::decompress::{ParDecompressBuilder}};

const INFO_FIELD_INDEX   : usize = 7;  /// 0-based expected column index of the INFO field.
const GENOTYPES_START_IDX: usize = 9;  /// 0-based expected column index where genotype entries are expected to begin.

/// Struct representing the VCF `INFO` field.
#[derive(Debug, Default)]
struct InfoField (Option<Vec<String>>);

impl Deref for InfoField {
    type Target = Vec<String>;
    fn deref(&self) -> &Self::Target {
        self.0.as_ref().expect("Attempting to access empty InfoField!")
    }
}

impl InfoField {
    /// Instantiate a new InfoField from a provided string slice.
    /// # Arguments:
    /// - `field`: raw, unparsed `INFO` vcf column string slice.
    pub fn new(field: &str) -> Self {
        Self(Some(field.split(';').map(|s| s.to_string()).collect::<Vec<String>>()))
    }

    /// Clear-out the structs content.
    pub fn clear(&mut self) {
        self.0 = None
    }

    /// return `true` if InfoField contains the "MULTI_ALLELIC" tag. (i.e. is multiallelic.) 
    pub fn is_multiallelic(&self) -> bool {
        self.iter().any(|field| field == "MULTI_ALLELIC")
    }

    /// return `true` if InfoField contains the "VT=SNP" tag. (i.e. is an SNP)
    pub fn is_snp(&self) -> Result<bool, String> {
        let vtype = self.iter()
            .find(|&field| field.starts_with("VT="))
            .ok_or_else(||{
                String::from("INFO field does not contain any 'VT=' tag.")
            })?
            .split('=')
            .collect::<Vec<&str>>()[1];
        Ok(vtype == "SNP")
    }

    /// Return the annotated population allele frequency for a given population.
    /// # Arguments:
    /// - `pop`: super-population-id. 
    /// 
    /// # Behavior: 
    /// for a given `pop` raw string slice, this method will search for any field matching "{pop}_AF="
    /// and return its value.
    pub fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, String> {
        let pop_af_regex = format!("{}_AF", pop);
        self.iter()
            .find(|&field| field.starts_with(&pop_af_regex))
            .map(|x| x.split('='))
            .ok_or_else(|| {
                return format!("INFO field does not contain any '{pop_af_regex}' tag.")
            })?
            .collect::<Vec<&str>>()[1]
            .parse::<f64>()
            .map_err(|err|{
                format!("Failed to parse population allele frequency to a float using '{pop_af_regex}'. Got ['{err}']")
            })
    }
}

/// GenotypeReader using a `.vcf`, of `.vcf.gz` file.
/// # Description 
/// Reader tries to be as efficient as possible, and will read lines as lazily as possible.  
/// Thus, lines are actually read field by field, as conservatively as possible, making various checks along the way
/// (e.g. 'is this chromosome coordinate relevant', 'is this line actually a bi-allelic SNP?', etc.) and skipped as soon as
/// a discrepancy is found out. 
/// 
/// # Fields:            
/// - `source`          : Boxed BufReader for the given `.vcf(.gz)` file.
/// - `samples`         : Vector of Sample-Id. These are extracted from the VCF Header, from fields 9 to n
/// - `info`            : `InfoField` representing the `INFO` column of the VCF file, for the current position.
/// - `buf`             : Raw string-byte buffer.
/// - `genotypes_filled`: boolean representing whether or not `buf` field is filled.
/// - `idx`             : index counter, indicating at which column index the reader is currently at.
/// 
/// # Implemented Traits: 
/// - `GenotypeReader`
pub struct VCFReader<'a> {
    pub source           : Box<dyn BufRead + 'a>,
    samples              : Vec<String>,
    info                 : InfoField,
    buf                  : Vec<u8>,
    pub genotypes_filled : bool,
    idx                  : usize,
}

impl<'a> VCFReader<'a> {
    /// Instantiate an initialize new VCFReader.
    /// # Arguments:
    /// - `path`: path leading to the `.vcf(.gz)` file.
    /// - `threads`: number of decompression threads (This is only relevant in the case of BGZF compressed `.vcf.gz` files)
    pub fn new(path: &Path, threads: usize) -> std::io::Result<VCFReader<'a>>{
        let mut reader = Self::get_reader(path, threads)?;
        let samples = Self::parse_samples_id(&mut reader)?;

        Ok(VCFReader{source: reader, samples, buf: Vec::new(), info: InfoField::default(), genotypes_filled: false, idx:0})
    }

    /// Skip to the next line field and fill the buffer with its contents.
    pub fn next_field(&mut self) -> Result<&str, Box<dyn Error>> {
        self.clear_buffer();
        self.source.read_until(b'\t', &mut self.buf)?;
        self.buf.pop();
        self.idx += 1;
        Ok(std::str::from_utf8(&self.buf)?)
    }

    /// Fill `self.buffer` until an EOL is encountered and reset the `self.idx` counter to 0.
    fn next_eol(&mut self) -> std::io::Result<()> {
        let _ = self.source.read_until(b'\n', &mut self.buf)?;
        self.idx=0;
        Ok(())
    }

    /// Skip to the next line and clear all buffers: `self.buf`, `self.info`, `self.idx`
    pub fn next_line(&mut self)-> std::io::Result<()> {
        if ! self.genotypes_filled {
            self.next_eol()?;
        }
        self.info.clear();
        self.clear_buffer();
        self.idx = 0; 
        self.genotypes_filled = false;
        Ok(())
    }

    /// Return `true` if the `INFO` field of the current line contains a `MULTI_ALLELIC` tag.
    pub fn is_multiallelic(&self) -> bool {
        self.info.is_multiallelic()
    }

    /// Return `true` if the `INFO` field of the current line contains a `VT=SNP` tag.
    pub fn is_snp(&self) -> Result<bool, String> {
        self.info.is_snp()
    }

    /// Skip `n` fields on the current line and increment `self.idx` accordingly.
    /// # Arguments: 
    /// - `n`: number of fields to skip. 
    pub fn skip(&mut self, n: usize) -> std::io::Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new())?;
        }
        self.idx+=n;
        Ok(())
    }

    /// Parse fields 0 and 1 of the VCF line and return the current chromosome and position coordinates.
    /// # Errors:
    /// - if `self.idx` != 0 (i.e. we're not currently at the beginning at the line.)
    pub fn parse_coordinate(&mut self) -> Result<(u8, u32), Box<dyn Error>> {
        if self.idx == 0 {
            let chromosome : u8  = self.next_field()?.parse().map_err(|err|{
                format!("'{err}' While attempting to parse chromosome coordinate.")
            })?; // 1
            let position   : u32 = self.next_field()?.parse().map_err(|err|{
                format!("'{err}' While attempting to parse position coordinate.")
            })?; // 2
            Ok((chromosome, position))
        } else {
            Err(format!("VCFReader line index is greater than 0: {}", self.idx).into())
        }
    }

    /// Skip to the expected `INFO` field idx and serialize it into `self.info` as a `InfoField` struct.
    /// # Errors:
    /// - if `self.idx > INFO_FIELD_INDEX constant.
    pub fn parse_info_field(&mut self) -> Result<(), Box<dyn Error>> {
        if self.idx <= INFO_FIELD_INDEX {
            self.skip(INFO_FIELD_INDEX - self.idx)?;
            self.info = InfoField::new(self.next_field()?);
            Ok(())
        } else {
            Err(format!("Cannot parse INFO field, as current field {} is greater than the expected INFO field index ({INFO_FIELD_INDEX})", self.idx).into())
        }
    }

    /// Skip to othe expected sample genotypes fields and fill `self.buf` with all subsequent fields until an EOL is found.
    /// - Once this method is called, `self.buf` is expected to contain all the samples genotypes for the current line.
    /// - Since we're working with raw bytes:
    ///   - REF (0) == 48
    ///   - ALT (1) == 49
    ///   - SEP (|) == 124
    pub fn fill_genotypes(&mut self) -> Result<(), Box<dyn Error>> {
        if self.idx <= GENOTYPES_START_IDX {
            self.clear_buffer();
            self.skip(GENOTYPES_START_IDX-self.idx)?;
            self.next_eol()?;
            self.genotypes_filled = true;
            Ok(())
        } else {
            Err(format!("Cannot parse genotype fields, as current field {} is greater than the expected genotypes fields start index ({GENOTYPES_START_IDX})", self.idx).into())
        }
    }

    /// Clear out the contents of `self.buf`
    pub fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    /// Returns `true` if the file has not yet encountered an EOF.
    pub fn has_data_left(&mut self) -> std::io::Result<bool> {
        self.source.fill_buf().map(|b| ! b.is_empty())
    }

    /// Clone and return the contents of `self.samples`
    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    /// Check the file extension of the provided file, and return an appropriate BufReader
    /// - `.vcf` -> Return a default BufReader
    /// - `.gz`  -> Return a parallel BGZF decompressor/reader
    /// # Arguments 
    /// - `path`   : path leading to the targeted vcf file.
    /// - `threads`: number of user-provided decompression threads for the BGZF decompressor.
    ///              (Only relevant if the file extension ends with `.gz`)
    fn get_reader(path: &Path, threads: usize) -> std::io::Result<Box<dyn BufRead>> {
        let path_ext = path.extension().ok_or_else(||{
            return std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid or missing file extension for vcf file: '{}'", path.display())
            )
        })?;
        let source: Result<Box<dyn Read>, std::io::Error> = match path_ext.to_str(){
            Some("vcf") => Ok(Box::new(File::open(path)?)),
            Some("gz")  => {
                let reader = File::open(path)?;
                let builder = ParDecompressBuilder::<Bgzf>::new()
                    .maybe_num_threads(threads)
                    .maybe_par_from_reader(reader);
                Ok(builder)
            },
            _           => {
                Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("Invalid file extension for vcf file: '{}'. \
                        Accepted format are ['vcf', 'vcf.gz']", path.display()
                    )
                ))
            }
        };
        Ok(Box::new(BufReader::new(source?)))
    }

    /// Skip all vcf description lines until the header line has been found (i.e. the line starts with '#CHROM').
    /// Then, fill `self.samples` with the contents of this line. 
    /// # Arguments:
    /// - `reader`: a BufReader targeting a vcf file.
    /// 
    /// # Panics:
    /// - If the reader reads all the file contents without encountering any line starting with the '#CHROM' pattern.
    fn parse_samples_id(reader: &mut Box<dyn BufRead + 'a>) -> std::io::Result<Vec<String>>{
        let mut samples = Vec::new();
        for line in reader.lines() {
            let line = line?;
            let split_line: Vec<&str> = line.split('\t').collect();
            if split_line[0] == "#CHROM" {
                for ind in &split_line[..]{
                    samples.push(ind.to_string());
                }
                return Ok(samples)
            }
        }
        panic!();
    }
}

impl<'a> GenotypeReader for VCFReader<'a> {
        // Return the alleles for a given SampleTag. Search is performed using `sample_tag.idx()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2], String> {
        let geno_idx = sample_tag.idx()
            .as_ref()
            .ok_or_else(|| {
                format!("Failed to obtain column index of sample {} while attempting to retrieve its alleles.", sample_tag.id())
            })? * 4;
            
        let retrieve_err = || format!("Failed to retrieve the alleles of sample {} within the current VCF.", sample_tag.id());
        let haplo1 = self.buf.get(geno_idx  ).ok_or_else(retrieve_err).map(|all| all - 48)?;
        let haplo2 = self.buf.get(geno_idx+2).ok_or_else(retrieve_err).map(|all| all - 48)?;
        Ok([haplo1, haplo2])
    }
    
    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, String> {
        Ok(self.info.get_pop_allele_frequency(pop)?)
    }
}