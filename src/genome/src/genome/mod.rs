use std::{
    collections::BTreeMap,
    ops::{Deref}, 
    path::Path,
    io::{BufRead, BufReader},
    fs::File, str::FromStr
};

use crate::{
    chromosome::Chromosome,
    coordinate::ChrIdx
};

mod error;
use error::FastaIndexReaderError;

use log::{warn, info, debug};

/// `BTreeMap` of chromosomes, with Key: chromosome name (u8) | Value: Chromosome
#[derive(Debug, Clone)]
pub struct Genome(BTreeMap<ChrIdx, Chromosome>);

impl Deref for Genome {
    type Target = BTreeMap<ChrIdx, Chromosome>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Genome {
    /// Instantiate a new, empty chromosome
    /// # @TODO convert this into default()
    #[must_use]
    pub fn new() -> Genome {
        Genome(BTreeMap::new())
    }

    /// Create a new genome from a slice of Chromosomes
    #[must_use]
    pub fn from(chromosomes: &[Chromosome]) -> Genome {
        let mut genome = Genome::new();
        for chr in chromosomes {
            genome.add_chromosome(*chr);
        }
        genome
    }

    /// Read a `.fasta.fai` file and parse it into a vector of Chromosome structs.
    /// # Arguments
    /// - `path`: Path leading to either a `.fasta`, `.fasta.fai`, `.fa` or `.fa.fai` file.
    ///           In the case of a `.fasta` or `.fa` file, a companion `.fai` with a matching file-stem
    ///           must be located within the same directory.
    /// 
    /// # Expected file format:
    /// - Fields         : `<CHROMOSOME>`    `<LENGTH>`
    /// - Field-separator: `'\t'`
    /// 
    /// # Errors
    /// - returns `FastaIndexReaderError::InvalidExt` if ...
    ///   - `path` does not carry a file extension.
    ///   - `path` extension is neither (`.fa`, `.fasta` or `.fai`)
    /// - returns `FastaIndexReaderError::FileNotFound` if no matching `.fai` file could be found within the 
    ///   target directory
    /// - returns `FastaIndexReaderError::ParseInt` if there is an invalid length value at column 1
    pub fn from_fasta_index(path: &str) -> Result<Genome,FastaIndexReaderError> {
        use FastaIndexReaderError::{InvalidExt, FileNotFound, ParseLine};
        info!("Parsing reference genome: {}", path);

        // ---- Extract the  file extention of `path` and format the expected `.fai` file if necessary.
        let Some(ext) = Path::new(path).extension().map(|ext| ext.to_string_lossy().to_string()) else {
            return Err(InvalidExt{ext:"None".to_string()})
        };
        let fai = match ext.as_str() { 
            "fai"          => path.to_string(),                          
            "fasta" | "fa" => format!("{}{}", path, ".fai"), // Append '.fai' if it is not there
            other          => return Err(InvalidExt{ext: other.to_owned()})
        };
        debug!("Matching fasta.fai file : {}", fai);
    
    
        let mut genome = Self::new();
        let file = BufReader::new(match File::open(fai.clone()) {
            Ok(file) => file,
            Err(err) => return Err(FileNotFound{path: fai, err: err.to_string()}),
        });
    
        let mut skipped_chrs = Vec::new();
        for (index,line) in file.lines().enumerate() {
            let line = line.map_err(|err| ParseLine{idx: index, err: err})?;
            match Chromosome::from_str(&line) {
                Err(_) => {
                    let skipped = line.split("\t").take(1).collect::<Vec<&str>>()[0].to_string();
                    skipped_chrs.push(skipped);
                }
                Ok(chr) =>{ 
                    genome.add_chromosome(chr);
                    debug!("Chromosome: {: <10} {: <12}", chr.name, chr.length);
                }            
            };
        }
    
        if !skipped_chrs.is_empty(){
            warn!("\nSome chromosomes were skipped while parsing {}:\n{:?}", fai, skipped_chrs);
        }
        Ok(genome)
    }

    pub fn add_chromosome(&mut self, chromosome: Chromosome) -> Option<Chromosome> {
        self.0.insert(chromosome.name, chromosome)
    }
}

/// Simply returns a default genome index in case the user did not provide a specific .fasta.fai file. 
impl Default for Genome {
    fn default() -> Self {
        warn!("No reference genome provided. Using default.");
        Genome(BTreeMap::from([
            (ChrIdx::from( 1), Chromosome::new( 1, 249_250_621)),
            (ChrIdx::from( 2), Chromosome::new( 2, 243_199_373)),
            (ChrIdx::from( 3), Chromosome::new( 3, 198_022_430)),
            (ChrIdx::from( 4), Chromosome::new( 4, 191_154_276)),
            (ChrIdx::from( 5), Chromosome::new( 5, 180_915_260)),
            (ChrIdx::from( 6), Chromosome::new( 6, 171_115_067)),
            (ChrIdx::from( 7), Chromosome::new( 7, 159_138_663)),
            (ChrIdx::from( 8), Chromosome::new( 8, 146_364_022)),
            (ChrIdx::from( 9), Chromosome::new( 9, 141_213_431)),
            (ChrIdx::from(10), Chromosome::new(10, 135_534_747)),
            (ChrIdx::from(11), Chromosome::new(11, 135_006_516)),
            (ChrIdx::from(12), Chromosome::new(12, 133_851_895)),
            (ChrIdx::from(13), Chromosome::new(13, 115_169_878)),
            (ChrIdx::from(14), Chromosome::new(14, 107_349_540)),
            (ChrIdx::from(15), Chromosome::new(15, 102_531_392)),
            (ChrIdx::from(16), Chromosome::new(16,  90_354_753)),
            (ChrIdx::from(17), Chromosome::new(17,  81_195_210)),
            (ChrIdx::from(18), Chromosome::new(18,  78_077_248)),
            (ChrIdx::from(19), Chromosome::new(19,  59_128_983)),
            (ChrIdx::from(20), Chromosome::new(20,  63_025_520)),
            (ChrIdx::from(21), Chromosome::new(21,  48_129_895)),
            (ChrIdx::from(22), Chromosome::new(22,  51_304_566))
        ]))
    }
}
