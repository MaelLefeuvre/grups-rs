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

const INFO_FIELD_INDEX   : usize = 7;
const GENOTYPES_START_IDX: usize = 9;

#[derive(Debug, Default)]
struct InfoField (Option<Vec<String>>);

impl Deref for InfoField {
    type Target = Vec<String>;

    fn deref(&self) -> &Self::Target {
        self.0.as_ref().expect("Attempting to access empty InfoField!")
    }
}

impl InfoField {
    pub fn new(field: &str) -> Self {
        Self(Some(field.split(';').map(|s| s.to_string()).collect::<Vec<String>>()))
    }

    pub fn clear(&mut self) {
        self.0 = None
    }

    pub fn is_multiallelic(&self) -> bool {
        self.iter().any(|field| field == "MULTI_ALLELIC")
    }

    pub fn is_snp(&self) -> bool {
        let vtype = self.iter()
            .find(|&field| field.starts_with("VT="))
            .unwrap()
            .split('=')
            .collect::<Vec<&str>>()[1];
        vtype == "SNP"
    }

    pub fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, std::num::ParseFloatError> {
        self.iter()
            .find(|&field| field.starts_with(&format!("{}_AF", pop)))
            .unwrap()
            .split('=')
            .collect::<Vec<&str>>()[1]
            .parse::<f64>()
    }
}

pub struct VCFReader<'a> {
    pub source       : Box<BufReader<Box<dyn Read + 'a>>>,
    samples          : Vec<String>,
    info             : InfoField,
    buf              : Vec<u8>,
    pub genotypes_filled : bool,
    idx              : usize,
}

impl<'a> VCFReader<'a> {
    pub fn new(path: &Path, threads: usize) -> std::io::Result<VCFReader<'a>>{
        let mut reader = Self::get_reader(path, threads)?;
        let samples = Self::parse_samples_id(&mut reader)?;

        Ok(VCFReader{source: reader, samples, buf: Vec::new(), info: InfoField::default(), genotypes_filled: false, idx:0})
    }

    pub fn next_field(&mut self) -> Result<&str, Box<dyn Error>> {
        self.clear_buffer();
        self.source.read_until(b'\t', &mut self.buf)?;
        self.buf.pop();
        self.idx += 1;
        Ok(std::str::from_utf8(&self.buf)?)
    }

    fn next_eol(&mut self) -> std::io::Result<()> {
        let _ = self.source.read_until(b'\n', &mut self.buf)?;
        self.idx=0;
        Ok(())
    }

    pub fn next_line(&mut self)-> std::io::Result<()> {
        self.info.clear();
        self.clear_buffer();
        self.idx = 0; 
        if ! self.genotypes_filled {
            self.next_eol()?;
        }
        self.genotypes_filled = false;
        Ok(())
    }

    pub fn is_multiallelic(&self) -> bool {
        self.info.is_multiallelic()
    }

    pub fn is_snp(&self) -> bool {
        self.info.is_snp()
    }

    pub fn skip(&mut self, n: usize) -> std::io::Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new())?;
        }
        self.idx+=n;
        Ok(())
    }

    pub fn parse_coordinate(&mut self) -> Result<(u8, u32), Box<dyn Error>> {
        if self.idx == 0 {
            let chromosome : u8  = self.next_field()?.parse().unwrap(); // 1
            let position   : u32 = self.next_field()?.parse().unwrap(); // 2
            Ok((chromosome, position))
        } else {
            Err(format!("VCFReader line index is greater than 0: {}", self.idx).into())
        }
    }

    pub fn parse_info_field(&mut self) -> Result<(), Box<dyn Error>> {
        if self.idx <= INFO_FIELD_INDEX {
            self.skip(INFO_FIELD_INDEX - self.idx)?;
            self.info = InfoField::new(self.next_field()?);
            Ok(())
        } else {
            Err(format!("Cannot parse INFO field, as current field {} is greater than the expected INFO field index ({INFO_FIELD_INDEX})", self.idx).into())
        }
    }

    pub fn fill_genotypes(&mut self) -> std::io ::Result<()> {
        self.clear_buffer();
        self.skip(GENOTYPES_START_IDX-self.idx)?;
        self.next_eol()?;
        self.genotypes_filled = true;
        Ok(())
    }

    pub fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    pub fn has_data_left(&mut self) -> std::io::Result<bool> {
        self.source.fill_buf().map(|b| ! b.is_empty())
    }

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    fn get_reader(path: &Path, threads: usize) -> std::io::Result<Box<BufReader<Box<dyn Read>>>> {
        let source: Box<dyn Read> = match path.extension().unwrap().to_str(){
            Some("vcf") => Box::new(File::open(path)?),
            Some("gz")  => {
                let reader = File::open(path)?;
                let builder = ParDecompressBuilder::<Bgzf>::new().maybe_num_threads(threads).maybe_par_from_reader(reader);
                Box::new(builder)
            }
            _           => panic!()
        };
        Ok(Box::new(BufReader::new(source)))
    }

    fn parse_samples_id(reader: &mut Box<BufReader<Box<dyn Read>>>) -> std::io::Result<Vec<String>>{
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
    fn get_alleles(&self, sample_tag: &SampleTag ) -> Option<[u8; 2]> {
        let geno_idx=sample_tag.idx().unwrap() * 4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        Some([haplo1, haplo2])
    }

    fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64, Box<dyn Error>> {
        Ok(self.info.get_pop_allele_frequency(pop)?)
    }
}