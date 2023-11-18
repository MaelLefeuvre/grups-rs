use std::{
    collections::{HashMap, BTreeMap},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    fs::File
};

use crate::{
    parse,
    read::{SampleTag, genotype_reader::VCFReader},
};

use log::warn;

use fastrand;
use anyhow::{Result, bail};
use located_error::{LocatedError, LocatedOption};

mod error;
pub use error::PanelReaderError;

pub const PANEL_EXT: [&str; 1] = ["panel"];

/// Iterate over the contents of a OS-directory and search for all files ending with `.panel`
/// Return these matches as a vector of PathBuf.
/// 
/// # Arguments:
/// - `input_dir`: path of the targeted directory, were search should be performed.
/// 
/// # Errors:
/// Returns an error if, after iterating over all the contents of `input_dir`:
/// - the output Vec<PathBuf> is empty.
/// - the length of the output Vec<PathBuf> is greater than 1. 
/// @TODO: use this as a 'from_dir()' constructor within PanelReader instead.
fn fetch_input_panel(input_dir: &Path) -> Result<PathBuf>{
    use PanelReaderError::{NotFound, MultipleFound};
    let panel = parse::fetch_input_files(input_dir, &PANEL_EXT)
        .loc("While Fetching input panel")?;

    let err_msg = || format!("While searching for a panel definition file in {}", input_dir.display());
    match panel.len() {
        1 => Ok(panel[0].clone()),
        0 => Err(NotFound).with_loc(err_msg),
        _ => Err(MultipleFound(panel)).with_loc(err_msg),
    }    
}




/// Input sample panel definition file reader (`.panel` extension)
/// ### File characteristics:
/// - Tab-separated fields
/// - Columns: <Sample-id>  <Super-pop-id>    <pop-id>    <gender (optional)>
/// ### Fields:
/// - source_file: path to the source `.panel` file
/// - samples: HashMap with K: (String) <pop-id> | V: Vec<SampleTag>
#[derive(Debug, Clone)]
pub struct PanelReader {
    pub source_file : PathBuf,
    pub samples: HashMap<String, Vec<SampleTag>>,
}


impl PanelReader {
    /// Instantiate a Panel-definition from a `.panel` file
    /// # Arguments:
    /// - `panel_path`: path leading to the `.panel` file.
    /// # TODO: Convert this into ::from()
    pub fn new(panel_path: &Path) -> Result<Self> {
        let source = BufReader::new(File::open(panel_path)?);
        let samples = Self::parse(source)?;
        Ok(Self{samples, source_file: panel_path.to_path_buf()})
    }    


    pub fn from_dir(target_dir: &Path) -> Result<Self> {
        let panel = fetch_input_panel(target_dir)?;
        Self::new(&panel)
    }

    /// Copy and subset the source `.panel` files.
    /// The output file will only contain entries that are contained within the `samples` field.
    /// # Fields
    /// - `output_dir`: path leading to the output `.panel` file.
    pub fn copy_from_source(&self, output_dir: &Path) -> Result<()> {
        let source     = BufReader::new(File::open(&self.source_file)?);
        let mut writer = BufWriter::new(File::create(output_dir)?);
        for line in source.lines(){
            let line = line?;
            let split_line = line.split('\t').collect::<Vec<&str>>();
            let (pop, superpop) = (split_line[1], split_line[2]);
            if self.samples.contains_key(superpop) || self.samples.contains_key(pop) {
                writer.write_all(format!("{line}\n").as_bytes())?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    /// Return a random subsample of contaminating individuals, from given populations and number of individuals.
    /// # Arguments:
    /// - `contam_pop`: Vector of (super-)population-ids. Each entry corresponds to a given pileup individual.
    /// - `contam_num_ind`: Vector of user-defined required-individuals. Each entry corresponds to a given pileup individual.
    /// 
    /// # Output:
    /// A vector of vectors of SampleTags.
    /// - Each vector[i] is of size contam_num_ind[i]
    /// - Each SampleTag within vector[i] has been randomly sampled (without replacement) from contam_pop[i]
    /// 
    /// # Behavior:
    /// - `When contam_pop.len()` and `contam_num_ind.len()` differ, indices and values are expected to wrap around,
    ///    according to the smallest vector length out of the two. 
    pub fn fetch_contaminants(&self, contam_pop: &[String], contam_num_ind: &Vec<usize>) -> Result<Vec<Vec<SampleTag>>> {
        let mut samples_contam_tags: Vec<Vec<SampleTag>> = Vec::new();

        // ---- Loop along `contam_num_ind` and populate `samples_contam_tags` with n sampletags.
        for num in contam_num_ind {
            let mut contam_ind_tags : Vec<SampleTag> = Vec::new();

            // ---- Wrap if contam_pop.len() < contam_num_ind.len() OR if contam_pop.len() > contam_num_ind.len()
            let contam_pop = &contam_pop[*num % std::cmp::min(contam_num_ind.len(), contam_pop.len())];

            // ---- Randomly sample n sampletags
            // ---- @TODO: This type of error handling is performed multiple times. stay DRY.
            for _ in 0..*num {
                let contam_sample_tag = self.random_sample(contam_pop, None)?
                    .with_loc(|| PanelReaderError::MissingContaminant(contam_pop.to_string()))?
                    .clone();

                contam_ind_tags.push(contam_sample_tag);
            }
            samples_contam_tags.push(contam_ind_tags);
        }
        Ok(samples_contam_tags)
    }

    /// Randomly sample a SampleTag, given a population-id. Optionally subsample our panel before sampling.
    /// # Arguments:
    /// - `pop`: (super-)Population-id string
    /// - `exclude`: Optional vector of sample tags to exclude from our sampling batch.
    ///    Thus, any provided sample tag cannot become a return value.
    pub fn random_sample(&self, pop: &str, exclude: Option<&Vec<&SampleTag>>) -> Result<Option<&SampleTag>> {
        // Extract the vector of candidate SampleTags using the provided population-id.
        // Bailout if `pop` does not match anything.
        let candidates = self.samples.get(pop).with_loc(|| PanelReaderError::MissingSample(pop.to_string()))?;

        // If there are sample tags to exclude, pre-filter out these SampleTags from `self.samples` before random sampling.
        let candidate = candidates.iter()
            .filter(|sample| !exclude.unwrap_or(&vec![]).contains(sample))
            .collect::<Vec<&SampleTag>>();

        // ---- Error out if there are no candidates left.
        if candidate.is_empty() { bail!(PanelReaderError::ExhaustedPanel) }
        
        // ---- Return a random sample
        let rng = fastrand::Rng::new();   
        Ok(candidate.get(rng.usize(0..candidate.len())).copied())
    }

    /// Subset `self.samples` with entries matching the provided `subset_pops`.
    /// # Arguments:
    /// - `subset_pops`: slice of (super-)population-id.
    pub fn subset_panel(&mut self, subset_pops: &[String]) {
        self.samples.retain(|key,_| subset_pops.contains(key));
    }

    /// Convert the `self.samples` HashMap<String, Vec<SampleTag>> into a transposed BTreeMap<SampleTag, String>.
    pub fn into_transposed_btreemap(&self) -> BTreeMap<&SampleTag, Vec<&str>> {
        let mut transposed_data: BTreeMap<&SampleTag, Vec<&str>> = BTreeMap::new();
        for (pop, samples) in self.samples.iter(){
            
            for sample in samples.iter() {
                transposed_data.entry(sample).or_default().push(pop.as_str())
                //transposed_data.insert(value.clone(), key.clone());
            }
        }
        transposed_data
    }

    #[inline]
    pub fn flat_values(&self) -> impl Iterator<Item = (& String, & SampleTag)> {
        self.samples.iter().flat_map(|(pop_tag, sample_tags)| {
            sample_tags.iter().map( move |sample_tag| {
                (pop_tag, sample_tag)
            })
        })
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (&String, &Vec<SampleTag>)> {
        self.samples.iter()
    }

    /// Return a list of unique population keys
    #[inline]
    pub fn pop_keys(&mut self) -> impl Iterator<Item = &String>{
        self.samples.keys()
    }

    /// Pretty self-explanatory: sort each Vec<SampleTag> found within `self.samples`
    pub fn sort_panel(&mut self){
        self.samples.iter_mut().for_each(|(_, tag_vec)| tag_vec.sort());
    }

    /// Parse the provided input `.panel` definition file. and populate the `self.samples`
    pub fn parse(source: BufReader<File>) -> std::io::Result<HashMap<String, Vec<SampleTag>>> {
        let mut output: HashMap<String, Vec<SampleTag>> = HashMap::new();
        for line in source.lines(){
            let line = line?;
            let line: Vec<&str> = line.split('\t').collect();
            let sample_idx = None;
            output.entry(line[1].into()).or_default().push(SampleTag::new(line[0], sample_idx)); // Key == Super-pop
            output.entry(line[2].into()).or_default().push(SampleTag::new(line[0], sample_idx)); // Key == Pop
        }
        Ok(output)
    }

    /// Read a VCF file, search through its header for samples indices, and assign a field-index to each SampleTag found within `self.samples`
    /// # Arguments:
    /// - `vcf`: Path to the  targeted `.vcf` file. 
    pub fn assign_vcf_indexes(&mut self, vcf: &Path)  -> Result<()> {
        let loc_msg = || format!("While attempting to assign vcf indexes from {}", vcf.display());
        let reader = VCFReader::new(vcf, 0).with_loc(loc_msg)?;
        let header = &reader.samples();

        // ---- The intersect between the input panel file and the vcf might not be perfect -> Keep a record of missing entries.
        let mut skipped_entries = Vec::new();

        // ---- Loop along lines and assign indices for each vcf entry.
        for sample_tags in self.samples.values_mut(){
            for sample_tag in sample_tags.iter_mut() {
                let sample_idx = match header.iter().position(|id| id == sample_tag.id()){
                    Some(idx) => idx - 9, // Correct for previous fields.
                    None => {skipped_entries.push(sample_tag.id()); continue}
                };
                sample_tag.set_idx(sample_idx);
            }
        }

        // ---- Warn the user if some samples were not found within the vcf file.
        if ! skipped_entries.is_empty() {
            skipped_entries.sort();
            skipped_entries.dedup();
            warn!("Some sample entries were not found in input vcfs:\n{:?}", skipped_entries);
        }
        
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn random_sampletag(){
        let mut samples = HashMap::new();
        let pop = String::from("EUR");
        samples.insert(pop.to_owned(), vec![
            SampleTag::new("HG00096", Some(0)),
            SampleTag::new("NA06984", Some(1)),
        ]);
        let exclude = vec![
            SampleTag::new("NA06984", Some(1)),
            SampleTag::new("HO0000", Some(2)),
            SampleTag::new("HO0001", Some(3))
        ];
        samples.insert(String::from("HO"), exclude.clone());

        let source_file = std::path::Path::new("/dev/null");
        let panel = PanelReader{samples, source_file: source_file.to_path_buf()};

        // Since "NA06984" is part of the excluded inds, it's expected he'll never be sampled.
        for _ in 0..1000 {
            let random = panel.random_sample(&pop, Some(&exclude.iter().collect()))
                .expect("Failed to obtain random sample")
                .expect("Missing random sample");
            assert_eq!(random.id(), "HG00096");
        }
    }
}