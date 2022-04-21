//use flate2::read::GzDecoder;
use std::io::{Read, BufRead, BufReader, Lines};
use std::error::Error;
use std::borrow::Borrow;
use std::path::{PathBuf, Path};
use std::fs::File;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use log::{warn, info, debug};
use rand::seq::SliceRandom;

use rust_htslib::bgzf;
use rust_htslib::tpool::ThreadPool;

#[derive(Debug)]
pub struct RecombinationRate {
    chr     : u8,
    pos     : u32,
    _rate    : f64,
    _map     : f64,
}


impl RecombinationRate {
    pub fn new(line: &[&str]) -> Result<RecombinationRate, Box<dyn Error>> {
        let chr = str::replace(line[0], "chr", "").parse::<u8>()?;
        let pos = line[1].parse::<u32>()?;
        let _rate = line[2].parse::<f64>()?;
        let _map = line[3].parse::<f64>()?;
        Ok(RecombinationRate{chr, pos, _rate, _map})
    }
}

impl PartialEq<RecombinationRate> for RecombinationRate {
    fn eq(&self, other: &Self) -> bool { 
        self.chr == other.chr && self.pos == other.pos
    }
}

impl Eq for RecombinationRate {}

impl Hash for RecombinationRate {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.pos.hash(state);
    }
}

impl Borrow<u32> for RecombinationRate {
    fn borrow(&self) -> &u32 {
        self.pos.borrow()
    }
}

#[derive(Debug)]
pub struct GenMapReader {
    inner: HashMap<u8,HashSet<RecombinationRate>>
}


impl GenMapReader {
    pub fn new(input_path: &PathBuf) -> Result<GenMapReader, Box<dyn Error>> {
        let map_paths = Self::fetch_genetic_maps(input_path)?;
        let mut inner = HashMap::new();
        for map in map_paths {
            let source = BufReader::new(File::open(map.as_path())?);
            Self::parse_genetic_map(source, &mut inner)?;
        }
        Ok(GenMapReader{inner})
    }

    pub fn get(&self, chr: &u8, pos: &u32) -> Option<&RecombinationRate> {
        self.inner[chr].get(pos)
    }

    fn parse_genetic_map(source: BufReader<File>, inner: &mut HashMap<u8, HashSet<RecombinationRate>>) -> Result<(), Box<dyn Error>>{
        let mut lines = source.lines();
        lines.next(); // Skip header. 
        for line in lines {
            let line = line?;
            let line = &line.split('\t').collect::<Vec<&str>>();
            let recomb_rate = RecombinationRate::new(line)?;
            inner.entry(recomb_rate.chr).or_insert_with(HashSet::new).insert(recomb_rate);
        }
        Ok(())
    }

    fn fetch_genetic_maps(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{

        let paths = std::fs::read_dir(input_dir)?;
        let maps = paths.filter_map(Result::ok)
            .filter_map(|d| d.path()
                .to_str()
                .and_then(|f|
                    if f.ends_with(".txt") { 
                        Some(d) 
                    } 
                    else { 
                        None
                    }
                )
                .map(|f| f.path())
            )
        .collect::<Vec<PathBuf>>();
        Ok(maps)
    }
}

#[derive(Debug, Clone)]
pub struct Sample {
    id: String,
    idx: usize,
}

impl Sample{
    pub fn new(id: &str, idx: usize) -> Sample {
        Sample{id: id.to_string(), idx}
    }

    pub fn id(&self) -> &String {
        &self.id
    }
    pub fn idx(&self) -> &usize {
        &self.idx
    }
}

impl PartialEq for Sample {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Sample {}

impl std::cmp::Ord for Sample {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.cmp(&other.id)
    }
}

impl std::cmp::PartialOrd for Sample {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
pub struct VCFPanelReader{
    pub samples: HashMap<String, Vec<Sample>>,
}

impl VCFPanelReader {
    pub fn new(panel_path: &Path, vcf: &Path, tpool: &ThreadPool) -> std::io::Result<VCFPanelReader> {
        let source = BufReader::new(File::open(panel_path)?);
        let header = VCFReader::new(vcf, tpool)?.samples();
        let samples = Self::parse_header(source, &header)?;
        Ok(VCFPanelReader{samples})
    }

    pub fn random_sample(&self, pop: &String) -> Option<&Sample> {
        self.samples[pop].choose(&mut rand::thread_rng())
    }

    pub fn parse_header(source: BufReader<File>, header: &[String]) -> std::io::Result<HashMap<String, Vec<Sample>>> {
        let mut output: HashMap<String, Vec<Sample>> = HashMap::new();

        for line in source.lines(){
            let line = line?;
            let line: Vec<&str> = line.split('\t').collect();
            let sample_idx = match header.iter().position(|id| id == line[0]){
                Some(idx) => idx,
                None => {warn!("Sample not found in input vcfs. Skipping:\n{:?}", line); continue}
            };
            output.entry(line[1].into()).or_insert(Vec::new()).push(Sample::new(line[0], sample_idx));
            output.entry(line[2].into()).or_insert(Vec::new()).push(Sample::new(line[0], sample_idx));

        }
        Ok(output)
    }
}

pub struct VCFReader<'a> {
    pub source: Box<dyn BufRead + 'a>,
    samples   : Vec<String>,
    buf       : Vec<u8>,
    idx       : usize,
}

impl<'a> VCFReader<'a> {
    pub fn new(path: &Path, tpool: &ThreadPool) -> std::io::Result<VCFReader<'a>>{
        let mut reader = Self::get_reader(path, tpool)?;
        let samples = Self::parse_samples(&mut reader)?;

        Ok(VCFReader{source: reader, samples, buf: Vec::new(), idx:0})
    }

    fn next_field(&mut self) -> Result<&str, Box<dyn Error>> {
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

    fn skip(&mut self, n: usize) -> std::io::Result<()> {
        for _ in 0..n {
            self.source.read_until(b'\t', &mut Vec::new())?;
        }
        self.idx+=n;
        Ok(())
    }
    fn fill_genotypes(&mut self) -> std::io::Result<()> {
        let genotypes_start_idx = 9;
        self.clear_buffer();
        self.skip(genotypes_start_idx-self.idx)?;
        self.next_eol()?;
        Ok(())
    }

    fn get_alleles(&mut self, idx: &usize) -> Result<(String, String), Box<dyn Error>> {
        let geno_idx=idx*4;
        let haplo1 = std::str::from_utf8(&[self.buf[geno_idx]])?.to_owned();
        let haplo2 = std::str::from_utf8(&[self.buf[geno_idx+2]])?.to_owned();
        Ok((haplo1, haplo2))
    }
    fn clear_buffer(&mut self) {
        self.buf.clear();
    }

    fn has_data_left(&mut self) -> std::io::Result<bool> {
        self.source.fill_buf().map(|b| !b.is_empty())
    }

    pub fn parse_sample(&mut self, idx: &usize) -> Result<(), Box<dyn Error>> {
        let mut i = 0;
        let mut genotypes: Vec<u8> = Vec::new();
        while self.has_data_left()? {
            let chr: u8  = self.next_field()?.parse()?;
            self.clear_buffer();
            let pos: u32 = self.next_field()?.parse()?;
            self.clear_buffer();
            self.fill_genotypes()?;
            let (haplo1, haplo2) = self.get_alleles(idx)?;
            if i % 500000 == 0 {
                println!("{} {: >2} {: >9} {} {}",
                    idx, chr, pos, haplo1, haplo2
                );

            }
            self.clear_buffer();
            genotypes.clear();

            i+=1;
        }
        Ok(())
    }

    pub fn lines(self) -> Lines<Box<dyn BufRead + 'a>> {
        self.source.lines()
    }

    pub fn samples(&self) -> Vec<String> {
        self.samples.clone()
    }

    fn get_reader(path: &Path, tpool: &ThreadPool) -> std::io::Result<Box<BufReader<Box<dyn Read>>>> {
        let mut source = bgzf::Reader::from_path(path).unwrap();
        source.set_thread_pool(tpool).unwrap();
        let source = Box::new(source);
        
        //let source: Box<dyn Read> = match path.extension().unwrap().to_str(){
            //Some("f") => Box::new(File::open(path)?),
            //Some("gz")  => Box::new(GzDecoder::new(File::open(path)?)),
            //Some("gz")  => Box::new(bgzf::Reader::from_path(path).unwrap()),
           // _           => panic!()
        //};
        Ok(Box::new(BufReader::new(source)))

    }
    fn parse_samples(reader: &mut Box<BufReader<Box<dyn Read>>>) -> std::io::Result<Vec<String>>{
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

use crate::pedigree::*;

pub fn pedigree_parser(path: &Path) -> std::io::Result<Pedigree> {
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
                pedigree.add_individual(&label, None)?;
            },
            Some(ParseMode::Relationships) => {
                let (offspring, parent1, parent2) = parse_pedline(line, "=repro(")?;
                pedigree.set_relationship(&offspring, (&parent1,&parent2))?;

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