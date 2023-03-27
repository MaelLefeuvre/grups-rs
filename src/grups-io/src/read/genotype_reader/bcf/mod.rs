use std::{path::{Path, PathBuf}};


use genome::coordinate::{Coordinate};
use located_error::{LocatedError};

use crate::{
    parse,
    read::SampleTag,
    read::genotype_reader::{GenotypeReader},
};

use anyhow::Result;
use log::debug;


const VCF_EXT: [&str; 2] = ["vcf", "vcf.gz"];


use rust_htslib::bcf::{record::{GenotypeAllele, Genotypes}, record::{Buffer, Info}, IndexedReader, header::HeaderView, Record, Records};
use rust_htslib::bcf::Read as BcfRead;

pub struct BcfReader<'a> {
    inner        : IndexedReader,
    pub source   : PathBuf,
    pub record   : Record,
    pub gen_buffer : Buffer,
    pub inf_buffer : Buffer,
    genotypes    : Option<Genotypes<'a, Buffer>>
}

impl<'a> GenotypeReader for BcfReader<'a> {
    fn get_alleles(&self, sample_tag: &SampleTag) -> Result<[u8;
2]> {
        let Some(sample_idx) = *sample_tag.idx() else {
            panic!()
        };
        let sample = self.genotypes.as_ref().unwrap().get(sample_idx);
        let alleles = match (sample[0], sample[1]) {
            (GenotypeAllele::Unphased(l), GenotypeAllele::Phased(r)) => [l as u8, r as u8],
            _ => panic!()
        };

        Ok(alleles)
    }

    fn get_pop_allele_frequency(&self,pop: &str) -> Result<f32> {
        let tag = format!("{pop}_AF");
        let af = self.record.info(tag.as_bytes()).float().unwrap().unwrap()[0];
        Ok(af)

    }

    fn fetch_input_files(input_dir: &Path) -> Result<Vec<PathBuf>>where Self:Sized {
        let vcfs = parse::fetch_input_files(input_dir, &VCF_EXT).loc("While searching for candidate vcf files.")?;
        debug!("Found the following vcf file candidates as input for the pedigree simulations: {:#?}", vcfs);
        Ok(vcfs)
    }
}

impl<'a> BcfReader<'a> {
    pub fn new(path: &Path, threads: usize) -> Result<BcfReader> {
        let mut inner = IndexedReader::from_path(path).unwrap();
        inner.set_threads(threads).unwrap();
        let  record = inner.empty_record();
        let genotypes = None;
        let gen_buffer = Buffer::new();
        let inf_buffer = Buffer::new();
        let source = path.to_path_buf();
        Ok(BcfReader{inner, source, record, gen_buffer, inf_buffer, genotypes})
    }

    /// Return an iterator over all the file's record
    pub fn records(&mut self) -> Records<IndexedReader>{
        self.inner.records()
    }

    /// Set the bcf to the next position, and fill the buffer.
    /// WARNING: 0-based coordinate.
    pub fn set_record(&mut self, coordinate: &Coordinate) {
        let chr = coordinate.chromosome.0.to_string().bytes().collect::<Vec<u8>>();
        let rid = self.inner.header().name2rid(&chr).unwrap();
        let start = (u32::from(coordinate.position) ) as u64;
        let end   = Some(start+1);
        self.inner.fetch(rid, start, end).unwrap();
        
        self.inner.read(&mut self.record);
        while !self.record.is_snp() { // There might be CNVs upstream...
            self.inner.read(&mut self.record);
        }
    }

    /// 
    pub fn fill_genotypes(&'a mut self) {
        let genotypes = self.record.genotypes().unwrap();
        self.genotypes = Some(genotypes)
    }

    /// Shared buffer genotypes.
    pub fn genotypes(&mut self) -> Genotypes<&mut Buffer> {
        self.record.genotypes_shared_buffer(&mut self.gen_buffer).unwrap()
    }

    /// Shared buffer Info
    pub fn info(&mut self, tag: &'a [u8]) -> Info<&mut Buffer> {
        self.record.info_shared_buffer(tag, &mut self.inf_buffer)
    }

    /// Return true if the current Record contains a VT=SNP tag.
    pub fn is_snp(&self) -> bool {
        let vt = self.record.info(b"VT").string().unwrap().unwrap();
        vt.iter().any(|vt|std::str::from_utf8(vt).unwrap() != "SNP")
    }

    /// Return true if 
    /// - the current Record contains a "MULTI_ALLELIC" flag
    /// - allele_count > 2
    pub fn is_multiallelic(&self) -> bool {
        self.record.info(b"MULTI_ALLELIC").flag().unwrap() || self.record.allele_count() > 2
    }
}

pub trait VariantFilter {
    fn is_snp(&self) -> bool;
    fn is_multiallelic(&self) -> bool;
}

impl VariantFilter for Record {
    /// Return true if the current Record contains a VT=SNP tag.
    fn is_snp(&self) -> bool {
        let vt = self.info(b"VT").string().unwrap().unwrap();
        vt.iter().any(|vt|std::str::from_utf8(vt).unwrap() != "SNP")
    }

    /// Return true if 
    /// - the current Record contains a "MULTI_ALLELIC" flag
    /// - allele_count > 2
    fn is_multiallelic(&self) -> bool {
        self.info(b"MULTI_ALLELIC").flag().unwrap() || self.allele_count() > 2
    }
}