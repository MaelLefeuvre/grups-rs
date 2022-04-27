use std::cell::RefMut;
//use flate2::read::GzDecoder;
use std::io::{Read, BufRead, BufReader, Lines};
use std::error::Error;
use std::path::{PathBuf, Path};
use std::fs::File;
use std::collections::{HashMap, HashSet};
use log::{warn, info, debug};
use pwd_from_stdin::genome::{SNPCoord, Genome};
use rand::seq::SliceRandom;

use rust_htslib::bgzf;
use noodles_bgzf::{AsyncReader};
use tokio::io::{self, AsyncBufRead};
use tokio::io::AsyncBufReadExt;


#[derive(Debug, Clone)]
pub struct SampleTag {
    id: String,
    idx: usize,
}

impl SampleTag {
    pub fn new(id: &str, idx: usize) -> SampleTag {
        SampleTag{id: id.to_string(), idx}
    }

    pub fn id(&self) -> &String {
        &self.id
    }
    pub fn idx(&self) -> &usize {
        &self.idx
    }
}

impl PartialEq for SampleTag {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for SampleTag {}

impl std::cmp::Ord for SampleTag {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.cmp(&other.id)
    }
}

impl std::cmp::PartialOrd for SampleTag {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub struct VCFPanelReader{
    pub samples: HashMap<String, Vec<SampleTag>>,
}

impl VCFPanelReader {
    pub async fn new(panel_path: &Path, vcf: &Path) -> std::io::Result<VCFPanelReader> {
        let source = BufReader::new(File::open(panel_path)?);
        let reader = VCFReader::new(vcf)?;
        let samples = Self::parse_header(source, &reader.samples())?;
        drop(reader);
        Ok(VCFPanelReader{samples})
    }

    pub fn random_sample(&self, pop: &String) -> Option<&SampleTag> {
        self.samples[pop].choose(&mut rand::thread_rng())
    }

    pub fn parse_header(source: BufReader<File>, header: &[String]) -> std::io::Result<HashMap<String, Vec<SampleTag>>> {
        let mut output: HashMap<String, Vec<SampleTag>> = HashMap::new();

        for line in source.lines(){
            let line = line?;
            let line: Vec<&str> = line.split('\t').collect();
            let sample_idx = match header.iter().position(|id| id == line[0]){
                Some(idx) => idx,
                None => {warn!("Sample not found in input vcfs. Skipping:\n{:?}", line); continue}
            };
            output.entry(line[1].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx));
            output.entry(line[2].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx));

        }
        Ok(output)
    }
}


pub struct VCFAsyncReader {
    pub source: AsyncReader<tokio::fs::File>,
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

    pub fn get_alleles(&mut self, idx: &usize) -> Result<(u8, u8), Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        Ok((haplo1, haplo2))
    }

    pub fn get_alleles2(&mut self, idx: &usize) -> Result<Option<[u8; 2]>, Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        let alleles = Some([haplo1, haplo2]);
        Ok(alleles)
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

    async fn get_asyncreader(path: &Path, threads: usize) -> std::io::Result<AsyncReader<tokio::fs::File>> {
        use tokio::fs::File;

        let builder = noodles_bgzf::AsyncReader::builder(File::open(path).await?).set_worker_count(threads);
        let source = builder.build();

        Ok(source)
    }

    async fn parse_samples_id(reader: &mut AsyncReader<tokio::fs::File>) -> std::io::Result<Vec<String>>{
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

pub struct VCFReader<'a> {
    pub source: Box<BufReader<Box<dyn Read + 'a>>>,
    samples   : Vec<String>,
    buf       : Vec<u8>,
    idx       : usize,
}

impl<'a> VCFReader<'a> {
    pub fn new(path: &Path) -> std::io::Result<VCFReader<'a>>{
        let mut reader = Self::get_reader(path)?;
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

    pub fn get_alleles(&mut self, idx: &usize) -> Result<(u8, u8), Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = self.buf[geno_idx]   - 48;
        let haplo2 = self.buf[geno_idx+2] - 48;
        Ok((haplo1, haplo2))
    }

    pub fn get_alleles2(&mut self, idx: &usize) -> Result<Option<[u8; 2]>, Box<dyn Error>> {
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



    //pub async fn lines(self) -> Lines<Box<dyn BufRead + 'a>> {
    //    self.source.lines().await?
    //}

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    fn get_reader(path: &Path) -> std::io::Result<Box<BufReader<Box<dyn Read>>>> {
        let source = noodles_bgzf::Reader::new(File::open(path)?);
        let source = Box::new(source);
        Ok(Box::new(BufReader::new(source)))


        //Ok(source);
        //let source = Box::new(AsyncBufRead::)
         

        //let source: Box<dyn Read> = match path.extension().unwrap().to_str(){
            //Some("f") => Box::new(File::open(path)?),<
            //Some("gz")  => Box::new(GzDecoder::new(File::open(path)?)),
            //Some("gz")  => Box::new(bgzf::Reader::from_path(path).unwrap()),
           // _           => panic!()
        //};
        //Ok(Box::new(BufReader::new(source)))

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
}

use crate::pedigree::*;

pub fn pedigree_parser<'a> (path: &'a Path, genome: &'a Genome) -> std::io::Result<Pedigree> {
    #[derive(Debug)]
    enum ParseMode {Individuals, Relationships, Comparisons}
    let mut parse_mode = None;
    let mut pedigree = Pedigree::new();
    let reader = BufReader::new(File::open(path)?);
    for line in reader.lines() {
        let line= line?;
        let line: Vec<&str> = line.split('#').collect();
        match line[0].chars().next() { // Skip comments and empty lines.
            Some('#') | None => continue,
            Some(_)          => (),
        }

        let mut change_parse_mode = |mode| {parse_mode = Some(mode); true};
        if let true = match line[0] {
            "INDIVIDUALS"   => change_parse_mode(ParseMode::Individuals),
            "RELATIONSHIPS" => change_parse_mode(ParseMode::Relationships),
            "COMPARISONS"   => change_parse_mode(ParseMode::Comparisons),
            _               => (false)
        }{continue}

        match parse_mode {
            Some(ParseMode::Individuals)   => {
                let label = line[0].to_string();
                pedigree.add_individual(&label, None, genome.clone())?;
            },
            Some(ParseMode::Relationships) => {
                let (offspring, parent1, parent2) = parse_pedline(line, "=repro(")?;
                pedigree.set_relationship(&offspring, (&parent1, &parent2))?;

            },
            Some(ParseMode::Comparisons)   => {
                let (label, ind1, ind2) = parse_pedline(line, "=compare(")?;
                pedigree.add_comparison(&label, (&ind1, &ind2))?;
            },
            None                           => continue
        };
    }
    Ok(pedigree)
}

fn parse_pedline (line: Vec<& str>, regex: &str) -> std::io::Result<(String, String, String)> {
    use std::io::ErrorKind::InvalidData;
    let mut temp=line[0].trim()
        .strip_suffix(')')
        .ok_or(InvalidData)?
        .split(regex)
        .map(|s| s.to_string());

    let ind = temp.next().ok_or(InvalidData)?;

    let parents: Vec<String> = temp.next()
        .ok_or(InvalidData)?
        .split(',')
        .map(|s| s.to_string())
        .collect();
    Ok((ind, parents[0].to_owned(), parents[1].to_owned()))
}

pub fn get_input_vcfs(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{
    let paths = std::fs::read_dir(input_dir)?;

    let mut vcfs = paths.filter_map(Result::ok)
        .filter_map(|d| d.path()
            .to_str()
            .and_then(|f|
                if f.ends_with(".vcf") || f.ends_with(".vcf.gz") { 
                    debug!("Found: {}", f);
                    Some(d) 
                } 
                else { 
                    debug!("Skipping: {}", f);
                    None
                }
            )
            .map(|f| f.path())
        )
    .collect::<Vec<PathBuf>>();
    vcfs.sort();
    info!("Found input vcf file candidates: {:#?}", vcfs);
    Ok(vcfs)
}

pub fn fetch_input_panel(input_dir: &PathBuf) -> std::io::Result<PathBuf>{

    let paths = std::fs::read_dir(input_dir)?;
    let panel = paths.filter_map(Result::ok)
        .filter_map(|d| d.path()
            .to_str()
            .and_then(|f|
                if f.ends_with(".panel") { 
                    Some(d) 
                } 
                else { 
                    None
                }
            )
            .map(|f| f.path())
        )
    .collect::<Vec<PathBuf>>();

    if panel.len() != 1 {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, 
            format!("Found multiple candidate Panel definition files: {:#?}\n
            Please specify the relevant file using '--panel'.
            Exiting.", panel)
        ))
    }

    info!("Found: {}", panel[0].to_str().unwrap_or("None"));
    Ok(panel[0].clone())
}