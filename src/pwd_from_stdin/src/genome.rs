
use crate::jackknife::*;

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::error::Error;
use std::io::{BufReader, BufRead};
use log::{info, warn, debug};
use std::ops::{Range, Deref, DerefMut};
use rust_lapper::{Interval, Lapper};

use rand::Rng;

/// A simple struct representing an SNPCoordinate position.
/// Note that :
///   - REF/ALT positions are optional, but are required if the user
///     requires known-sites filtration.
///   - SNPCoordinates implement Equality with JackknifeBlock structs.
///     See: `pwd_from_stdin::jackknife::JackknifeBlock`
///   - Hashable, but only in regards to chr and pos.
///     -> don't insert alternate variations.
#[derive(Debug, Clone)]
pub struct SNPCoord {
    pub chromosome : u8,
    pub position   : u32,
    pub reference  : Option<char>,
    pub alternate  : Option<char>,
}

impl Ord for SNPCoord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.chromosome, self.position).cmp(&(other.chromosome, other.position))
    }
}

impl PartialOrd for SNPCoord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq<SNPCoord> for SNPCoord {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.position == other.position
    }
}

impl PartialEq<JackknifeBlock> for SNPCoord {
    fn eq(&self, other: &JackknifeBlock) -> bool {
        other.chromosome == self.chromosome && self.position >= other.range.start && self.position < other.range.end
    }
}

impl Eq for SNPCoord {}

impl Hash for SNPCoord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.position.hash(state);
    }
}

pub struct Allele {
    pos: u32,
    allele: u8,
    af: f64,
}

#[derive(Default)]
pub struct Alleles (Vec<Allele>);

impl Alleles{
    pub fn add_allele(&mut self, allele: Allele){
        self.0.push(allele);
    }
}

impl Extend<Allele> for Alleles {
    fn extend<T: IntoIterator<Item=Allele>>(&mut self, iter: T) {
        for elem in iter {
            self.add_allele(elem);
        }
    }
}

#[derive(Debug, Clone)]
pub struct Locus {
    pos: u32,
    alleles: (u8, u8),
    af: f64,
}

impl Locus {
    pub fn crossover(&mut self) {
        self.alleles = (self.alleles.1, self.alleles.0)

    }

    pub fn alleles(&self) -> (u8, u8) {
        self.alleles
    }
}

impl PartialEq<Locus> for Locus {
    fn eq(&self, other: &Self) -> bool { 
        self.pos == other.pos && self.alleles == other.alleles
    }
}

impl Eq for Locus {}

impl PartialOrd for Locus {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Locus {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.pos).cmp(&(other.pos))
    }
}

pub struct Chromatid {
    pub index   : usize, 
    pub name    : u8,
    pub length  : u32,
    pub alleles : Alleles
}

impl Chromatid {
    pub fn new(index: usize, name: u8, length: u32, alleles: Alleles) -> Self {
        Chromatid{index, name, length, alleles}
    }
}

#[derive(Debug, Clone)]
/// A simple struct representing a chromosome index. This is mainly used to compute Jackknife Blocks
pub struct Chromosome {
    pub index  : usize,
    pub name   : u8,
    pub length : u32,
    loci       : Vec<Locus>,

}

impl Chromosome {
    pub fn new(index: usize, name: u8, length: u32) -> Chromosome{
        Chromosome{index, name, length, loci: Vec::new()}
    }

    pub fn snp_len(&self) -> usize {
        self.loci.len()
    }

    pub fn loci(&self) -> impl Iterator<Item = &Locus> {
        self.loci.iter()
    }

    pub fn add_locus(&mut self, pos: u32, alleles: (u8, u8), af: f64) {
        self.loci.push(Locus{pos, alleles, af})
    }

    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }

    pub fn meiosis(&self, genetic_map: &GeneticMap) -> Chromatid {
        println!("Performing meosis in chr {}", self.name);

        let mut rng = rand::thread_rng();

        let mut previous_position = 0;
        let mut currently_recombining = false;

        let mut gamete = self.loci.clone();
        
        for locus in gamete.iter_mut() {
            let mut interval_prob_recomb: f64 = 0.0;

            for recombination_range in genetic_map[&self.name].find(previous_position, locus.pos) {
                let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
                let real_stop  = if locus.pos         > recombination_range.stop  {recombination_range.stop } else {locus.pos        };

                interval_prob_recomb += recombination_range.val.prob * (real_stop as f64 - real_start as f64 + 1.0);
            }

            if rng.gen::<f64>() < interval_prob_recomb {
                println!("Switch!");
                currently_recombining = ! currently_recombining;
            }

            if currently_recombining {
                locus.crossover();
            }

            previous_position = locus.pos;
        }

        let alleles: (Alleles, Alleles) = gamete.into_iter()
            .map(|Locus{pos, alleles, af}| {
                (Allele{pos, allele: alleles.0, af}, Allele{pos, allele: alleles.1, af})
            })
            .unzip();

        let random_chromatid = if rng.gen::<f64>() < 0.5 {alleles.0} else {alleles.1};
        Chromatid::new(self.index, self.name, self.length, random_chromatid)
    }
}

impl PartialEq for Chromosome {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}
impl Eq for Chromosome {}

impl PartialOrd for Chromosome {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chromosome {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.name).cmp(&other.name)
    }
}

impl Extend<Locus> for Chromosome {
    fn extend<T: IntoIterator<Item=Locus>>(&mut self, alleles: T) {
        for locus in alleles {
            self.add_locus(locus.pos, locus.alleles, locus.af);
        }
    }
}

#[derive(Debug)]
pub enum FastaIndexReaderError {
    FileNotFoundError(String, String),
    ParseIntError(String, u8, String),
    InvalidFileFormatError(String)
}
impl Error for FastaIndexReaderError {}

impl std::fmt::Display for FastaIndexReaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::FileNotFoundError(path, err)   => write!(f, "{}: {}",path, err),
            Self::ParseIntError(path, line, err) => write!(f, "{}: Line {} - Invalid line format [{}]", path, line, err),
            Self::InvalidFileFormatError(ext)  => write!(f, "Invalid fasta file format: expected [.fasta|.fa], got '{}'", ext),
        }
    }
}


#[derive(Default)]
pub struct Gamete(BTreeMap<u8, Chromatid>);

impl Deref for Gamete {
    type Target = BTreeMap<u8, Chromatid>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Gamete {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Gamete {

    pub fn add_chromatid(&mut self, chr: u8, chromatid: Chromatid) -> bool {
        self.insert(chr,  chromatid).is_none()
    }

    pub fn fertilize(&self, other: &Self) -> Genome {
        let mut genome = Genome::new();
        for ((_, chromatid1), (_, chromatid2)) in self.iter().zip(other.iter()) {
            assert_eq!(chromatid1.index, chromatid2.index);
            assert_eq!(chromatid1.name, chromatid2.name);
            assert_eq!(chromatid1.length, chromatid2.length);
                        assert_eq!(chromatid1.length, chromatid2.length);


            genome.add_chromosome(chromatid1.index, chromatid1.name, chromatid1.length);

            //let mut chromosome = Chromosome::new(chromatid1.index, chromatid1.name, chromatid1.length);
            for (allele1, allele2) in chromatid1.alleles.0.iter().zip(chromatid2.alleles.0.iter()) {
                assert_eq!(allele1.pos, allele2.pos);
                assert_eq!(allele1.af,  allele2.af );

                genome.get_chr_mut(&chromatid1.name).unwrap().add_locus(allele1.pos, (allele1.allele, allele2.allele), allele1.af)
            }
        }
        genome
    }
}

#[derive(Debug, Clone)]
pub struct Genome(BTreeMap<u8, Chromosome>);

impl Deref for Genome {
    type Target = BTreeMap<u8, Chromosome>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Genome {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Genome {

    fn new() -> Genome {
        Genome(BTreeMap::new())
    }

    

    pub fn from(chromosomes: &[Chromosome]) -> Genome {
        let mut genome = Genome(BTreeMap::new());
        for chr in chromosomes {
            genome.add_chromosome(chr.index, chr.name, chr.length);
        }
        genome
    }

    /// Read a `.fasta.fai` file and parse it into a vector of Chromosome structs.
    /// 
    /// TODO: - convert 'path' variable to &str
    ///       - Check for .fai.fai extension doubles.
    pub fn from_fasta_index(path: &str) -> Result<Genome,FastaIndexReaderError> {
        use FastaIndexReaderError::{InvalidFileFormatError, FileNotFoundError, ParseIntError};
        info!("Parsing reference genome: {}", path);

        let fai = match Path::new(path).extension() { // Check file extension
            None                 => return Err(InvalidFileFormatError("None".to_string())),
            Some(osstr)          => match osstr.to_str(){ 
                Some("fai")          => path.to_string(),                          
                Some("fasta" | "fa") => format!("{}{}", path, ".fai"), // Append '.fai' if it is not there
                Some(other)    => return Err(InvalidFileFormatError(".".to_string()+other)),
                None                 => return Err(InvalidFileFormatError("None".to_string())),
            }        
        };
        debug!("fasta.fai file : {}", fai);
    
    
        let mut genome = Self::new();
        let file = BufReader::new(match fs::File::open(fai.clone()) {
            Ok(file) => file,
            Err(err) => return Err(FileNotFoundError(fai, err.to_string())),
        });
    
        let mut skipped_chrs = Vec::new();
        for (index,line) in file.lines().enumerate() {
            let line = line.unwrap();
            let split_line: Vec<&str> = line.split('\t').collect();
            match split_line[0].replace("chr", "").parse::<u8>() {
                Ok(result) =>{ 
                    let name      : u8        = result;
                    let length    : u32       = match split_line[1].parse() {
                        Ok(len) => len,
                        Err(e)  => return Err(ParseIntError(fai, name, e.to_string()))
                    };
                    debug!("Chromosome: {: <3} {: <10} {: <12}", index, name, length);
                    genome.add_chromosome(index, name, length);
                }
                Err(_) => skipped_chrs.push(split_line[0].to_string())
            };
        }
    
        if !skipped_chrs.is_empty(){
            warn!("\nSome chromosomes were skipped while parsing {}: {:#?}", fai, skipped_chrs);
        }
        Ok(genome)
    }

    fn add_chromosome(&mut self, index: usize, name: u8, length: u32) -> bool {
        self.0.insert(name, Chromosome::new(index, name, length)).is_none()
    }

    pub fn get_chr_mut(&mut self, name: &u8) -> Option<&mut Chromosome> {
        self.0.get_mut(name)
    }

    pub fn is_empty(&self) -> bool {
        self.0.values().all(|chr| chr.is_empty())
    }


    pub fn meiosis(&mut self, genetic_map: &GeneticMap) -> Gamete {
        let mut gamete = Gamete::default();
        for chr in self.0.values_mut() {
            if ! chr.is_empty(){
                let chromatid = chr.meiosis(genetic_map);
                gamete.add_chromatid(chr.name, chromatid);
            }
        }
        gamete
    }

}

/// Simply returns a default genome index in case the user did not provide a specific .fasta.fai file. 
impl Default for Genome {
    fn default() -> Self {
        warn!("No reference genome provided. Using default.");
        Genome(BTreeMap::from([
            (1, Chromosome::new( 0,  1, 249250621)),
            (2, Chromosome::new( 1,  2, 243199373)),
            (3, Chromosome::new( 2,  3, 198022430)),
            (4, Chromosome::new( 3,  4, 191154276)),
            (5, Chromosome::new( 4,  5, 180915260)),
            (6, Chromosome::new( 5,  6, 171115067)),
            (7, Chromosome::new( 6,  7, 159138663)),
            (8, Chromosome::new( 7,  8, 146364022)),
            (9, Chromosome::new( 8,  9, 141213431)),
            (10, Chromosome::new( 9, 10, 135534747)),
            (11, Chromosome::new(10, 11, 135006516)),
            (12, Chromosome::new(11, 12, 133851895)),
            (13, Chromosome::new(12, 13, 115169878)),
            (14, Chromosome::new(13, 14, 107349540)),
            (15, Chromosome::new(14, 15, 102531392)),
            (16, Chromosome::new(15, 16,  90354753)),
            (17, Chromosome::new(16, 17,  81195210)),
            (18, Chromosome::new(17, 18,  78077248)),
            (19, Chromosome::new(18, 19,  59128983)),
            (20, Chromosome::new(19, 20,  63025520)),
            (21, Chromosome::new(20, 21,  48129895)),
            (22, Chromosome::new(21, 22,  51304566))
        ]))
    }
}

#[derive(Debug, Clone)]
pub struct RecombinationRange {
    range: Range<u32>,
    prob : f64,
}

impl RecombinationRange {
    pub fn new(start: u32, end: u32, rate: f64) -> RecombinationRange {
        let prob: f64 = rate / 100.0 / 1_000_000.0 ;  // rate/cM/Mb.
        RecombinationRange{range: Range{start, end}, prob}
    }

    pub fn prob(&self) -> &f64 {
        &self.prob
    }

    pub fn rate(&self) -> f64 {
        self.prob * 100.0 * 1_000_000.0
    }
}

impl PartialEq<RecombinationRange> for RecombinationRange {
    fn eq(&self, other: &Self) -> bool { 
        self.range == other.range
    }
}

impl PartialEq<Range<u32>> for RecombinationRange {
    fn eq(&self, other: &Range<u32>) -> bool {
        self.range == *other
    }
}

impl Eq for RecombinationRange {}

impl Hash for RecombinationRange {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.range.hash(state);
    }
}

impl std::borrow::Borrow<Range<u32>> for RecombinationRange {
    fn borrow(&self) -> &Range<u32> {
        self.range.borrow()
    }
}

#[derive(Default)]
pub struct GeneticMap(HashMap<u8, Lapper<u32, RecombinationRange>>);

impl Deref for GeneticMap {
    type Target = HashMap<u8, Lapper<u32, RecombinationRange>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for GeneticMap {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl GeneticMap {
    pub fn from_dir(mut self, dir: &PathBuf) -> Result<GeneticMap, Box<dyn Error>> {
        let map_paths = Self::fetch_genetic_maps(dir)?;
        for map in map_paths.iter() {
            self.from_map(map)?;
        }
        Ok(self)
    }
    pub fn from_map(&mut self, path: &Path) -> Result<(), Box<dyn Error>> {
        let source = BufReader::new(File::open(path)?);
        let mut intervals: HashMap<u8, Vec<Interval<u32, RecombinationRange>>> = HashMap::new();

        let mut start = 0;
        for line in source.lines().skip(1) { // Skip header
            let line= line?;
            let line = &line.split('\t').collect::<Vec<&str>>();

            let chr   = str::replace(line[0], "chr", "").parse::<u8>()?;
            let stop = line[1].parse::<u32>()?;
            let rate = line[2].parse::<f64>()?;

            let recomb_rate = RecombinationRange::new(start, stop, rate);
            let interval = Interval{start, stop, val: recomb_rate};
            intervals.entry(chr).or_insert_with(Vec::new).push(interval);
            start = stop;
        }
        for (chr, intervals) in intervals.into_iter(){
            self.insert(chr, Lapper::new(intervals));
        }
        Ok(())
    }

    fn fetch_genetic_maps(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{

        let paths = std::fs::read_dir(input_dir)?;
        let maps = paths.filter_map(Result::ok)
            .filter_map(|d| d.path()
                .to_str()
                .and_then(|f|
                    if f.ends_with(".txt") {Some(d)} else {None}
                )
                .map(|f| f.path())
            )
        .collect::<Vec<PathBuf>>();
        Ok(maps)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_full_equality() {
        let chromosome1 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        let chromosome2 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }

    #[test]
    fn snp_partial_equality() {
        let chromosome1 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('T'), alternate: Some('G')};
        let chromosome2 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }
}
