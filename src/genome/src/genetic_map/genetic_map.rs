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

    pub fn compute_recombination_prob(&self, chromosome: u8, previous_position: u32, current_position: u32) -> f64 {
        let mut interval_prob_recomb = 0.0;
        for recombination_range in self[&chromosome].find(previous_position, current_position) {
            let real_start = if previous_position < recombination_range.start {recombination_range.start} else {previous_position};
            let real_stop  = if current_position  > recombination_range.stop  {recombination_range.stop } else {current_position };

            interval_prob_recomb += recombination_range.val.prob() * (real_stop as f64 - real_start as f64 + 1.0);
        }

        interval_prob_recomb
    }
}

#[cfg(test)]
mod tests {

}
