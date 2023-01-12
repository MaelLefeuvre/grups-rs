use std::{collections::HashMap, path::{Path, PathBuf}, io::{BufReader, BufRead}, fs::{self, File}};

use rust_lapper::{Interval, Lapper};

use located_error::prelude::*;
//use anyhow::{anyhow, Result};

mod error;
use error::*;

mod recomb_range;
use recomb_range::RecombinationRange;

use crate::coordinate::{Coordinate, Position, ChrIdx};


/// `HashMap` of Genetic Maps, in the form of a BITS Tree (Key = chromosome name | Value = Interval Tree)
/// 
/// See: <https://doi.org/10.1093/bioinformatics/bts652>
#[derive(Default)]
pub struct GeneticMap(HashMap<ChrIdx, Lapper<u32, RecombinationRange>>);

impl GeneticMap {
    /// Instantiate a `GeneticMap` from an OS directory.
    /// 
    /// # Arguments
    /// - `dir`: path leading to a directory containing genetic recombination maps (`.txt`)
    /// 
    /// # Errors
    /// - if no genetic map was found within the input directory
    pub fn from_dir(dir: impl AsRef<Path>) -> Result<Self> {
        let mut out = Self::default();
        for map in Self::fetch_genetic_maps(dir)? {
            out.from_map(&map).with_loc(|| GeneticMapError::ParseMap{map: map.clone()})?;
        }

        match out.0.keys().len() {
            0 => Err(anyhow::anyhow!(GeneticMapError::EmptyDir)),
            _ => Ok(out)
            // 22 => Ok(out)
        }.loc(GeneticMapError::EmptyDir)
    }

    /// Parse a BITS tree from a given genetic recombination map file, and add it to `self.0`
    /// Expected fields of a recombination map file are:
    ///  <Chromosome>    <Position(bp)>    <Recomb. rate (cM/Mb)>    <Map(cM)>
    /// 
    /// # Arguments
    /// - `path`: path leading to a genetic recombination map
    /// 
    /// # Errors
    /// - when the value of `path` is not found and/or does not have read permissions
    /// - can return either `ParseIntError` or `ParseFloatError` if one of the fields
    ///   contains invalid information. 
    pub fn from_map(&mut self, path: impl AsRef<Path> ) -> Result<()> {
        use GeneticMapError::*;

        // ---- Open the genetic recombination file and initialize HashMap
        let source = BufReader::new(File::open(path).loc("Failed to open file")?);
        let mut intervals: HashMap<ChrIdx, Vec<Interval<u32, RecombinationRange>>> = HashMap::new();

        // ---- Parse file and generate intervals from each line.
        let mut start = 0;
        for (i, line) in source.lines().enumerate().skip(1) { // Skip header
            let line = line.with_loc(|| InvalidLine(i))?;
            if line.is_empty(){ continue } // Skip empty lines.
            let fields: Vec<&str> = line.split('\t').collect::<Vec<&str>>();

            match fields.len() { // baild if we don't have 4 fields.
                4 => Ok(()),
                n => Err(anyhow!("Expected 4 fields, got {n}")).with_loc(|| InvalidFields(i))
            }?;

            let chr : ChrIdx = fields[0].parse().with_loc(|| ParseChr(i))?;
            let stop: u32    = fields[1].parse::<u32>().with_loc(|| ParsePos(i))?;
            let rate: f64    = fields[2].parse::<f64>().with_loc(|| ParseRate(i))?;

            // ---- Instantiate a new Interval, tied to a recombination range.
            let interval = Interval{ start, stop, val: RecombinationRange::new(start, stop, rate) };
            intervals.entry(chr).or_default().push(interval);
            start = stop;
        }

        // ---- Internalize intervals within `self`
        self.0.extend(intervals.into_iter().map(|(c,int)| (c, Lapper::new(int))));

        Ok(())
    }

    /// Iterate through the contents of an OS directory and return a list of found `.txt` files.
    /// # Arguments
    /// - `input_dir`: path leading to a directory containing genetic recombination maps (`.txt`)
    fn fetch_genetic_maps(input_dir: impl AsRef<Path>) -> Result<impl Iterator <Item = PathBuf>> {
        // Get a list of files within the specified directory
        let paths = fs::read_dir(input_dir).with_loc(|| GeneticMapError::ReadDir)?
            .filter_map(Result::ok)        
            .filter(|f| {
                // Filter out anything that does not end with the .txt file extension (case insensitive).
                f.path().extension().map_or(false, |ext| ext.eq_ignore_ascii_case("txt"))
            }).map(|f| f.path());
        
        Ok(paths)
    }

    /// Compute a probability of genetic recombination on a given chromosome, within a range
    /// (typically, the current and the previous position)
    /// # Parameters
    /// - `chromosome`        : name of the chromosome where the probability should be computed 
    /// - `previous_positions`: 0-based coordinate of the previously typed position.
    /// - `current_position`  : 0-based coordinate of the current position.
    /// 
    /// # Panics:
    /// - if `chromosome` does not match any key within `self`
    #[must_use]
    pub fn compute_recombination_prob(&self, coordinate: &Coordinate, previous_position: Position) -> f64 {
        // @TODO: Lapper Struct only allows for items implementing num_traits::PrimInt and num_trait::Unsigned.
        //This conversion could be avoided if I manage to implement these traits for position.
        let current_position : u32 = coordinate.position.into(); 
        let previous_position: u32 = previous_position.into();
        let mut interval_prob_recomb = 0.0;
        // ---- Search for all intervals contained between the range [previous_position, current_position[
        for recombination_range in self.0[&coordinate.chromosome].find(previous_position, current_position) {
            let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
            let real_stop  = if current_position  > recombination_range.stop  {recombination_range.stop } else {current_position };

            interval_prob_recomb += recombination_range.val.prob() * (f64::from(real_stop) - f64::from(real_start) + 1.0);
        }

        // Probability that an odd number of cross-over occurs, under Poisson distribution.
        0.5 * (1.0 - f64::exp(-2.0 * interval_prob_recomb))
    }
}