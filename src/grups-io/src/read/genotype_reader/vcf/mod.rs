use std::{io::{self, BufRead, BufReader, Read}, path::{Path, PathBuf}, fs::File};

mod info;
use info::InfoField;

use genome::coordinate::{Coordinate, Position, ChrIdx};
use located_error::{LocatedError, LocatedOption};
use located_error::loc;

use crate::{
    parse,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader, GenotypeReaderError},
};

use gzp::{deflate::Bgzf, par::decompress::ParDecompressBuilder};
use anyhow::Result;
use log::debug;

mod error;
pub use error::VCFReaderError;

const INFO_FIELD_INDEX   : usize = 7;  /// 0-based expected column index of the INFO field.
const GENOTYPES_START_IDX: usize = 9;  /// 0-based expected column index where genotype entries are expected to begin.
const VCF_EXT: [&str; 2] = ["vcf", "vcf.gz"];


impl<'a> GenotypeReader for VCFReader<'a> {
    // Return the alleles for a given SampleTag. Search is performed using `sample_tag.idx()`;
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Result<[u8; 2]> {
        use GenotypeReaderError::{MissingAlleles, InvalidSampleIndex};

        let sample_genotypes: Vec<&[u8]> = self.buf.split(|c| *c == b'\t' || *c == b'\n').collect();

        let geno_idx = sample_tag.idx().expect("Missing sample tag index"); 
        let sample_genotype = sample_genotypes.get(geno_idx).ok_or(InvalidSampleIndex)?;

        //trace!("  - {} | {}:{}", Coordinate::try_from(&self.coordinate_buffer[0..5])?, sample_tag.id(), std::str::from_utf8(sample_genotype)?);

        let retrieve_err = || format!("Failed to retrieve the alleles of sample {} within the current VCF.", sample_tag.id());
        let (haplo1, haplo2) = match sample_genotype.len() {
            3 => Ok(( // Autosomal, Pseudo-autosomal region or female X-chromosome
                sample_genotype.first().ok_or(MissingAlleles).map(|all| all - 48).with_loc(retrieve_err)?,
                sample_genotype.get(2).ok_or(MissingAlleles).map(|all| all - 48).with_loc(retrieve_err)?
            )),

            1 => { // Male X-chromosome
                let haplo = sample_genotype.first().ok_or(MissingAlleles).map(|all| all - 48).with_loc(retrieve_err)?;
                Ok((haplo, haplo))
            }
            _ => Err(MissingAlleles)
        }?;

        //let geno_idx = 4 * sample_tag.idx().as_ref().ok_or(InvalidSampleIndex)
        //    .with_loc(|| format!("While retrieving alleles of {}", sample_tag.id()))?;
            
        //let haplo1 = left //self.buf.get(geno_idx  ).ok_or(MissingAlleles).map(|all| all - 48).with_loc(retrieve_err)?;
        //let haplo2 = right //self.buf.get(geno_idx+2).ok_or(MissingAlleles).map(|all| all - 48).with_loc(retrieve_err)?;
        Ok([haplo1, haplo2])
    }
    
    // Return the alleles frequencies for a given population id.
    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f32> {
        use GenotypeReaderError::MissingFreq;
        self.info.get_pop_allele_frequency(pop)
            .map_err(|_|MissingFreq(pop.to_string()))
            .loc("While parsing coordinate")
    }

    fn fetch_input_files(input_dir: &Path) -> Result<Vec<PathBuf>> {
        let vcfs = parse::fetch_input_files(input_dir, &VCF_EXT).loc("While searching for candidate vcf files.")?;
        debug!("Found the following vcf file candidates as input for the pedigree simulations: {:#?}", vcfs);
        Ok(vcfs)
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
    pub fn new(path: &'a Path, threads: usize) -> Result<VCFReader<'a>>{
        let loc_msg = "While attempting to create a new VCFReader";
        let mut reader = Self::get_reader(path, threads).loc(loc_msg)?;
        let samples = Self::parse_samples_id(&mut reader).loc(loc_msg)?;

        Ok(VCFReader{source: reader, samples, buf: Vec::new(), info: InfoField::default(), genotypes_filled: false, idx:0})
    }

    /// Skip to the next line field and fill the buffer with its contents.
    pub fn next_field(&mut self) -> Result<&str> {
        use VCFReaderError::{FillBuffer, InvalidField};
        let loc_msg = "While attempting to parse next field" ;
        self.clear_buffer();
        self.source.read_until(b'\t', &mut self.buf).map_err(FillBuffer).loc(loc_msg)?;
        self.buf.pop();
        self.idx += 1;
        str::from_utf8(&self.buf).map_err(InvalidField).loc(loc_msg)
    }

    /// Fill `self.buffer` until an EOL is encountered and reset the `self.idx` counter to 0.
    fn next_eol(&mut self) -> Result<()> {
        let _ = self.source.read_until(b'\n', &mut self.buf)
            .map_err(VCFReaderError::FillBuffer)
            .loc("While attempting to find next EOL")?;
        self.idx=0;
        Ok(())
    }

    /// Skip to the next line and clear all buffers: `self.buf`, `self.info`, `self.idx`
    pub fn next_line(&mut self)-> Result<()> {
        if ! self.genotypes_filled {
            self.next_eol().map_err(VCFReaderError::SkipLineError)?;
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
    pub fn is_snp(&self) -> Result<bool> {
        self.info.is_snp().loc("While checking if this coordinate is an SNP")
    }

    /// Skip `n` fields on the current line and increment `self.idx` accordingly.
    /// # Arguments: 
    /// - `n`: number of fields to skip. 
    pub fn skip(&mut self, n: usize) -> Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new())
                .map_err(VCFReaderError::FillBuffer)
                .with_loc(|| format!("While attempting to skip {n} fields within the line"))?;
        }
        self.idx+=n;
        Ok(())
    }

    /// Parse fields 0 and 1 of the VCF line and return the current chromosome and position coordinates.
    /// # Errors:
    /// - if `self.idx` != 0 (i.e. we're not currently at the beginning at the line.)
    pub fn parse_coordinate(&mut self) -> Result<Coordinate> {
        use VCFReaderError::{ParseChrError, ParsePosError, InvalidFieldIdx};
        let loc_msg = "While parsing the genomic coordinate of this entry";
        if self.idx != 0 {
            return Err(InvalidFieldIdx(0, self.idx)).loc(loc_msg)
        }
        let chromosome : ChrIdx   = self.next_field()?.parse().map_err(ParseChrError).loc(loc_msg)?; // 1
        let position   : Position = self.next_field()?.parse().map_err(ParsePosError).loc(loc_msg)?; // 2
        Ok(Coordinate::new(chromosome, position))
    }

    /// Skip to the expected `INFO` field idx and serialize it into `self.info` as a `InfoField` struct.
    /// # Errors:
    /// - if `self.idx > INFO_FIELD_INDEX constant.
    pub fn parse_info_field(&mut self) -> Result<()> {
        let loc_msg = "While attempting to parse INFO field";
        if self.idx > INFO_FIELD_INDEX {
            return Err(VCFReaderError::InvalidFieldIdx(INFO_FIELD_INDEX, self.idx)).loc(loc_msg)
        }
        
        self.skip(INFO_FIELD_INDEX - self.idx).loc(loc_msg)?;
        self.info = InfoField::new(self.next_field().loc(loc_msg)?);
        Ok(())
    }

    /// Skip to the expected sample genotypes fields and fill `self.buf` with all subsequent fields until an EOL is found.
    /// - Once this method is called, `self.buf` is expected to contain all the samples genotypes for the current line.
    /// - Since we're working with raw bytes:
    ///   - REF (0) == 48
    ///   - ALT (1) == 49
    ///   - SEP (|) == 124
    pub fn fill_genotypes(&mut self) -> Result<()> {
        let loc_msg = "While attempting to parse INFO field";
        if self.idx > GENOTYPES_START_IDX {
            return Err(VCFReaderError::InvalidFieldIdx(GENOTYPES_START_IDX, self.idx)).loc(loc_msg)
        }
        self.clear_buffer();
        self.skip(GENOTYPES_START_IDX-self.idx).loc(loc_msg)?;
        self.next_eol().loc(loc_msg)?;
        self.genotypes_filled = true;
        Ok(())
    }

    /// Clear out the contents of `self.buf`
    pub fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    /// Returns `true` if the file has not yet encountered an EOF.
    pub fn has_data_left(&mut self) -> io::Result<bool> {
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
    ///   (Only relevant if the file extension ends with `.gz`)
    fn get_reader(path: &Path, threads: usize) -> Result<Box<dyn BufRead>> {
        use VCFReaderError::{InvalidFileExt, Open};
        let path_ext = path.extension().with_loc(||InvalidFileExt)?;
        let vcf      = File::open(path).with_loc(|| Open)?;
        let source: Box<dyn Read> = match path_ext.to_str(){
            Some("vcf") => Box::new(vcf),
            Some("gz")  => ParDecompressBuilder::<Bgzf>::new().maybe_num_threads(threads).maybe_par_from_reader(vcf),
            _           => return loc!(InvalidFileExt)
        };
        Ok(Box::new(BufReader::new(source)))
    }

    /// Skip all vcf description lines until the header line has been found (i.e. the line starts with '#CHROM').
    /// Then, fill `self.samples` with the contents of this line. 
    /// # Arguments:
    /// - `reader`: a BufReader targeting a vcf file.
    /// 
    /// # Panics:
    /// - If the reader reads all the file contents without encountering any line starting with the '#CHROM' pattern.
    fn parse_samples_id(reader: &mut Box<dyn BufRead + 'a>) -> Result<Vec<String>>{
        let mut samples = Vec::new();
        for line in reader.lines() {
            let line = line.with_loc(|| VCFReaderError::ParseLine)?;
            let split_line: Vec<&str> = line.split('\t').collect();
            if split_line[0] == "#CHROM" {
                samples.reserve(split_line.len());
                samples.extend(split_line.iter().map(|ind| ind.to_string()));
                return Ok(samples)
            }
        }
        panic!("Missing VCF-header within file.");
    }
}

#[cfg(test)]
mod tests {
    use gzp::{par::compress::{ParCompress, ParCompressBuilder}, ZWriter};

    use super::*;
    use std::io::Write;
    const FAKE_VCF: &str = "\
    ##fileformat=VCFv4.1\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\n\
    12\t60020\t.\tT\tTA,TAC\t100\tPASS\tAC=10,92;AMR_AF=0.0029,0.0086;AFR_AF=0.0008,0.0635;EUR_AF=0,0.002;SAS_AF=0.0031,0;EAS_AF=0.004,0;VT=INDEL;MULTI_ALLELIC\tGT\t0|0\t0|0\t0|0\t0|0\n\
    12\t60026\t.\tA\tC\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.5;AFR_AF=0.5;EUR_AF=0.5;SAS_AF=0.5;EAS_AF=0.5;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\n\
    12\t60057\t.\tC\tA\t100\tPASS\tAC=1582;AF=0.31;AMR_AF=0.02;AFR_AF=0.04;EUR_AF=0.06;SAS_AF=0.08;EAS_AF=0.010;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\n\
    12\t60083\t.\tG\tA\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.43;AFR_AF=0.38;EUR_AF=0.99;SAS_AF=0.61;EAS_AF=0.01;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1|1\n\
    ";

    const FAKE_VCF_XCHR: &str = "\
    ##fileformat=VCFv4.1\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\n\
    X\t60020\t.\tT\tTA,TAC\t100\tPASS\tAC=10,92;AMR_AF=0.0029,0.0086;AFR_AF=0.0008,0.0635;EUR_AF=0,0.002;SAS_AF=0.0031,0;EAS_AF=0.004,0;VT=INDEL;MULTI_ALLELIC\tGT\t0\t0|0\t0|0\t0|0\n\
    X\t60026\t.\tA\tC\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.5;AFR_AF=0.5;EUR_AF=0.5;SAS_AF=0.5;EAS_AF=0.5;AA=.|||;VT=SNP\tGT\t0|0\t0|1\t1|0\t1\n\
    X\t60057\t.\tC\tA\t100\tPASS\tAC=1582;AF=0.31;AMR_AF=0.02;AFR_AF=0.04;EUR_AF=0.06;SAS_AF=0.08;EAS_AF=0.010;AA=.|||;VT=SNP\tGT\t0\t0|1\t1|0\t1|1\n\
    X\t60083\t.\tG\tA\t100\tPASS\tAC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.43;AFR_AF=0.38;EUR_AF=0.99;SAS_AF=0.61;EAS_AF=0.01;AA=.|||;VT=SNP\tGT\t0\t0|1\t1\t1|1\n\
    ";

    #[test]
    fn test_open_vcf() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf");
        let mut file = File::create(&vcf_path)?;
        writeln!(file, "{FAKE_VCF}")?;

        let reader = VCFReader::new(&vcf_path, 1);
        assert!(reader.is_ok());
        Ok(())
    }

    #[test]
    fn test_open_vcf_gz() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf.gz");

        
        let file = File::create(&vcf_path)?;
        let mut parz: ParCompress<Bgzf> = ParCompressBuilder::new().from_writer(file);
        parz.write_all(FAKE_VCF.as_bytes()).expect("Failed to write VCF with ParCompressBuilder");
        parz.finish().expect("ParCompress should be able to flush its output.");

        let reader = VCFReader::new(&vcf_path, 1);
        assert!(reader.is_ok());
        Ok(())
    }

    #[test]
    fn test_open_vcf_invalid() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("README.txt");
        let _ = File::create(&vcf_path)?;
        let reader = VCFReader::new(&vcf_path, 1);
        assert!(reader.is_err());
        Ok(())
    }

    #[test]
    fn test_parse_sample_id() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf");
        let mut file = File::create(&vcf_path)?;
        writeln!(file, "{FAKE_VCF}")?;

        let reader = VCFReader::new(&vcf_path, 1).expect("Failed to create test reader");

        let expected_samples=["HG00096", "HG00097", "HG00099", "HG00100"];
        assert_eq!(reader.samples()[9..], expected_samples);
        Ok(())
    }

    #[test]
    fn test_readlines_autosomes() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf");
        let mut file = File::create(&vcf_path)?;
        writeln!(file, "{FAKE_VCF}")?;

        let mut reader = VCFReader::new(&vcf_path, 1).expect("Failed to create test reader");

        assert_eq!(reader.parse_coordinate()?, Coordinate::new(12, 60020));
        reader.parse_info_field()?;
        assert!(reader.is_multiallelic());
        reader.fill_genotypes()?;
        assert!(reader.genotypes_filled);

        reader.next_line()?;
        assert_eq!(reader.parse_coordinate()?, Coordinate::new(12, 60026));

        reader.next_line()?;
        assert_eq!(reader.parse_coordinate()?, Coordinate::new(12, 60057));

        reader.next_line()?;
        assert_eq!(reader.parse_coordinate()?, Coordinate::new(12, 60083));

        Ok(())
    }

    #[test]
    fn test_get_alleles_autosomes() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf");
        let mut file = File::create(&vcf_path)?;
        writeln!(file, "{FAKE_VCF}")?;

        let mut reader = VCFReader::new(&vcf_path, 1).expect("Failed to create test reader");

        let samples = &reader.samples()[9..];

        let expected_alleles = [
            [[0u8, 0], [0, 0], [0, 0], [0, 0]], // GT   0|0 0|0 0|0 0|0
            [[0, 0], [0, 1], [1, 0], [1, 1]],   // GT   0|0 0|1 1|0 1|1
            [[0, 0], [0, 1], [1, 0], [1, 1]],   // GT   0|0 0|1 1|0 1|1
            [[0, 0], [0, 1], [1, 0], [1, 1]]    // GT   0|0 0|1 1|0 1|1
        ];
        for allele_row in expected_alleles {
            assert_eq!(reader.parse_coordinate()?.chromosome, ChrIdx(12));

            reader.fill_genotypes()?;
            for (want, sample) in allele_row.iter().zip(samples.iter().enumerate().map(|(i, s)| SampleTag::new(s, Some(i), None))) {
                let got = reader.get_alleles(&sample)?;
                println!("{sample} ({want:?}) {got:?}");
                assert_eq!(*want, got)
            }
            reader.next_line()?;
        }
        Ok(())
    }

    #[test]
    fn test_get_alleles_x_chr() -> Result<()> {
        let tmpdir = tempfile::tempdir()?;
        let vcf_path = tmpdir.path().join("panel.vcf");
        let mut file = File::create(&vcf_path)?;
        writeln!(file, "{FAKE_VCF_XCHR}")?;

        let mut reader = VCFReader::new(&vcf_path, 1).expect("Failed to create test reader");

        let samples = &reader.samples()[9..];


        let expected_alleles = [
            [[0u8, 0], [0, 0], [0, 0], [0, 0]], // GT   0   0|0 0|0 0|0
            [[0, 0], [0, 1], [1, 0], [1, 1]],   // GT   0|0 0|1 1|0 1
            [[0, 0], [0, 1], [1, 0], [1, 1]],   // GT   0   0|1 1|0 1|1
            [[0, 0], [0, 1], [1, 1], [1, 1]]    // GT   0   0|1 1   1|1
        ];
        for allele_row in expected_alleles {
            assert_eq!(reader.parse_coordinate()?.chromosome, ChrIdx(b'X'));

            reader.fill_genotypes()?;
            for (want, sample) in allele_row.iter().zip(samples.iter().enumerate().map(|(i, s)| SampleTag::new(s, Some(i), None))) {
                let got = reader.get_alleles(&sample)?;
                println!("{sample} ({want:?}) {got:?}");
                assert_eq!(*want, got)
            }
            reader.next_line()?;
        }
        Ok(())
    }
}