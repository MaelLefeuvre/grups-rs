use std::{
    path::Path,
    pin::Pin,
    error::Error

};

use tokio::io::{AsyncBufRead, AsyncBufReadExt};


pub struct VCFAsyncReader {
    pub source: Pin<Box<dyn AsyncBufRead>>,
    samples   : Vec<String>,
    buf       : Vec<u8>,
    idx       : usize,
}

impl VCFAsyncReader {
    pub async fn new(path: &Path, threads: usize) -> std::io::Result<VCFAsyncReader>{
        let mut reader = Self::get_asyncreader(path, threads).await?;
        let samples = Self::parse_samples_id(&mut reader).await?;

        Ok(VCFAsyncReader{source: reader, samples, buf: Vec::new(), idx:0})
    }

    pub async fn next_field(&mut self) -> Result<&str, Box<dyn Error>> {
        self.clear_buffer();
        self.source.read_until(b'\t', &mut self.buf).await?;
        self.buf.pop();
        self.idx += 1;
        Ok(std::str::from_utf8(&self.buf)?)
    }

    async fn next_eol(&mut self) -> std::io::Result<()> {
        let _ = self.source.read_until(b'\n', &mut self.buf).await?;
        self.idx=0;
        Ok(())
    }

    pub async fn skip_line(&mut self) -> std::io::Result<()>{
        self.next_eol().await?;
        self.clear_buffer();
        Ok(())
    }


    pub async fn skip(&mut self, n: usize) -> std::io::Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new()).await?;
        }
        self.idx+=n;
        Ok(())
    }
    pub async fn fill_genotypes(&mut self) -> std::io ::Result<()> {
        let genotypes_start_idx = 9;
        self.clear_buffer();
        self.skip(genotypes_start_idx-self.idx).await?;
        self.next_eol().await?;
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
    


    pub fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    pub async fn has_data_left(&mut self) -> std::io::Result<bool> {
        Ok(! self.source.fill_buf().await?.is_empty())
    
    }

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    async fn get_asyncreader(path: &Path, threads: usize) -> std::io::Result<Pin<Box<dyn tokio::io::AsyncBufRead>>> {
        use tokio::fs::File;
        let source: Pin<Box<dyn tokio::io::AsyncBufRead>> = match path.extension().unwrap().to_str(){
            Some("vcf") => Box::pin(Box::new(tokio::io::BufReader::new(File::open(path).await?))),
            Some("gz")  => Box::pin(Box::new(noodles_bgzf::AsyncReader::builder(File::open(path).await?).set_worker_count(threads).build())),
            _           => panic!()
        };
        Ok(source)
    }

    async fn parse_samples_id(reader: &mut Pin<Box<dyn tokio::io::AsyncBufRead>>) -> std::io::Result<Vec<String>>{
        let mut samples = Vec::new();
        let mut lines = reader.lines();
        while let Some(line) = lines.next_line().await.unwrap() {
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