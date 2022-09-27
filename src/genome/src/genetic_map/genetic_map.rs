use std::{
    collections::{HashMap},
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
    error::Error,
    io::{BufReader, BufRead},
    fs::File,
};

use rust_lapper::{Interval, Lapper};

use super::recomb_range::RecombinationRange;

/// `HashMap` of Genetic Maps, in the form of a BITS Tree (Key = chromosome name | Value = Interval Tree)
/// 
/// See: <https://doi.org/10.1093/bioinformatics/bts652>
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
    /// Instantiate a `GeneticMap` from an OS directory.
    /// 
    /// # Arguments
    /// - `dir`: path leading to a directory containing genetic recombination maps (`.txt`)
    /// 
    /// # Errors
    /// - if no genetic map was found within the input directory
    pub fn from_dir(mut self, dir: &PathBuf) -> Result<GeneticMap, Box<dyn Error>> {
        let map_paths = Self::fetch_genetic_maps(dir)?;
        if map_paths.is_empty() {
            let dir_str = dir.to_str().expect("Failed to parse genetic-map directory to string.");
            return Err(format!("Failed to find or parse any genetic-map in provided directory: {dir_str}", ).into())
        }
        for map in &map_paths {
            self.from_map(map).map_err(|err| {
                format!(" Got [{err}] When attempting to parse genetic map '{map:?}'. \
                    Ensure the provided '--recomb-dir' ONLY contains properly formatted recombination map files."
                )
            })?;
        }
        Ok(self)
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
    pub fn from_map(&mut self, path: &Path) -> Result<(), Box<dyn Error>> {
        // ---- Open the genetic recombination file and initialize HashMap
        let source = BufReader::new(File::open(path)?);
        let mut intervals: HashMap<u8, Vec<Interval<u32, RecombinationRange>>> = HashMap::new();

        // ---- Parse file and generate intervals from each line.
        let mut start = 0;
        for line in source.lines().skip(1) { // Skip header
            let line= line?;
            if line.is_empty(){ continue } // Skip empty lines.
            let line = &line.split('\t').collect::<Vec<&str>>();
            let chr   = str::replace(line[0], "chr", "").parse::<u8>()?;
            let stop = line[1].parse::<u32>()?;
            let rate = line[2].parse::<f64>()?;

            // ---- Instantial a new Interval, tied to a recombination range.
            let recomb_rate = RecombinationRange::new(start, stop, rate);
            let interval = Interval{start, stop, val: recomb_rate};
            intervals.entry(chr).or_insert_with(Vec::new).push(interval);
            start = stop;
        }

        // ---- Internalize intervals within `self`
        for (chr, intervals) in intervals {
            self.insert(chr, Lapper::new(intervals));
        }
        Ok(())
    }

    /// Iterate through the contents of an OS directory and return a list of found `.txt` files.
    /// # Arguments
    /// - `input_dir`: path leading to a directory containing genetic recombination maps (`.txt`)
    fn fetch_genetic_maps(input_dir: &PathBuf) -> std::io::Result<Vec<PathBuf>>{

        // Get a list of files within the specified directory
        let paths = std::fs::read_dir(input_dir)?;

        // Filter out anything that does not end with the .txt file extension (case insensitive).
        let maps = paths.filter_map(Result::ok)
            .filter(|d| {
                d.path().extension()
                .map_or(false, |ext|  ext.eq_ignore_ascii_case("txt"))
            })
            .map(|f| f.path())
            .collect::<Vec<PathBuf>>();

        Ok(maps)
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
    pub fn compute_recombination_prob(&self, chromosome: u8, previous_position: u32, current_position: u32) -> f64 {
        let mut interval_prob_recomb = 0.0;
        // ---- Search for all intervals contained between the range [previous_position, current_position[
        for recombination_range in self[&chromosome].find(previous_position, current_position) {
            let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
            let real_stop  = if current_position  > recombination_range.stop  {recombination_range.stop } else {current_position };

            interval_prob_recomb += recombination_range.val.prob() * (f64::from(real_stop) - f64::from(real_start) + 1.0);
        }

        // Probability that an odd number of cross-over occurs, under Poisson distribution.
        0.5 * (1.0 - f64::exp(-2.0 * interval_prob_recomb))
    }
}
