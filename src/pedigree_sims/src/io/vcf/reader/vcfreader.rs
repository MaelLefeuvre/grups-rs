use std::{
    io::{BufRead, BufReader, Read},
    path::Path,
    fs::File,
    error::Error, cell::RefMut,
    collections::HashSet,
};

use genome::SNPCoord;

use crate::pedigree::Individual;

use gzp::{deflate::Bgzf, par::decompress::{ParDecompressBuilder}};


pub struct VCFReader<'a> {
    pub source: Box<BufReader<Box<dyn Read + 'a>>>,
    samples   : Vec<String>,
    buf       : Vec<u8>,
    idx       : usize,
}

impl<'a> VCFReader<'a> {
    pub fn new(path: &Path, threads: usize) -> std::io::Result<VCFReader<'a>>{
        let mut reader = Self::get_reader(path, threads)?;
        let samples = Self::parse_samples_id(&mut reader)?;

        Ok(VCFReader{source: reader, samples, buf: Vec::new(), idx:0})
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

    pub fn skip_line(&mut self) -> std::io::Result<()>{
        self.next_eol()?;
        self.clear_buffer();
        Ok(())
    }


    pub fn skip(&mut self, n: usize) -> std::io::Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new())?;
        }
        self.idx+=n;
        Ok(())
    }
    pub fn fill_genotypes(&mut self) -> std::io ::Result<()> {
        let genotypes_start_idx = 9;
        self.clear_buffer();
        self.skip(genotypes_start_idx-self.idx)?;
        self.next_eol()?;
        Ok(())
    }

    pub fn get_alleles(&self, idx: &usize) -> Result<(u8, u8), Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        Ok((haplo1, haplo2))
    }

    pub fn get_alleles2(&self, idx: &usize) -> Result<Option<[u8; 2]>, Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        let alleles = Some([haplo1, haplo2]);
        Ok(alleles)
    }


    pub fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    pub fn has_data_left(&mut self) -> std::io::Result<bool> {
        Ok(self.source.fill_buf().map(|b| ! b.is_empty())?)
    
    }

    pub fn parse_samples(&mut self, mut samples: Vec<RefMut<'_, Individual>>, valid_positions: &HashSet<SNPCoord>, pop: &str) -> Result<(), Box<dyn Error>> {
        let mut i = 0;
        while self.has_data_left()? {


            let chromosome : u8  = self.next_field()?.parse()?; // 1
            let position   : u32 = self.next_field()?.parse()?; // 2
            
            if i % 50_000 == 0 {
                println!("{i: >9} {chromosome: >2} {position: >9}");
            }
            i+=1;

            if ! valid_positions.contains(&SNPCoord{chromosome, position, reference: None, alternate: None}){
                self.skip_line()?;
                continue
            }

            self.skip(5)?;                                      // 6 
            let info = self.next_field()?.split(';').collect::<Vec<&str>>();

            if info.iter().any(|&field| field == "MULTI_ALLELIC") {
                self.skip_line()?;
                continue
            }
            
            let vtype = info.iter()
                .find(|&&field| field.starts_with("VT=")).unwrap()
                .split('=')
                .collect::<Vec<&str>>()[1];

            let pop_af = match vtype {
                "SNP" => {
                    let af = info.iter()
                    .find(|&&field| field.starts_with(&format!("{pop}_AF")))
                    .unwrap()
                    .split('=')
                    .collect::<Vec<&str>>()[1]
                    .parse::<f64>().unwrap();
                    Some(af)
                },
                _ => None
            };

            self.fill_genotypes()?;
            for founder in samples.iter_mut() {
                let alleles = self.get_alleles(founder.get_tag().unwrap().idx())?;
                founder.add_locus(&chromosome, position, alleles, pop_af.unwrap())?;
            }
        }

        for sample in samples{
            format!("{}: {}", sample.label, sample.genome[&1].snp_len());
        }

        Ok(())
    }

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    fn get_reader(path: &Path, threads: usize) -> std::io::Result<Box<BufReader<Box<dyn Read>>>> {


        let source: Box<dyn Read> = match path.extension().unwrap().to_str(){
            Some("vcf") => Box::new(File::open(path)?),
            //Some("gz")  => Box::new(noodles_bgzf::Reader::new(File::open(path)?)),
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
        let mut lines = reader.lines();
        while let Some(line) = lines.next() {
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

    pub fn compute_local_cont_af(&self, contam_ind_ids: &Vec<usize>) -> Result<f64, Box<dyn Error>>{
        let mut ref_allele_count = 0;
        let mut alt_allele_count = 0;
        for idx in contam_ind_ids.iter() {
            let cont_alleles = self.get_alleles2(idx)?.unwrap();
            for allele in cont_alleles.into_iter() {
                match allele {
                    0 => ref_allele_count +=1,
                    1 => alt_allele_count +=1,
                    _ => panic!("Contaminating individual is multiallelic.")
                }
            }
        }
        Ok((alt_allele_count as f64) /(alt_allele_count as f64 + ref_allele_count as f64))
    }

}
