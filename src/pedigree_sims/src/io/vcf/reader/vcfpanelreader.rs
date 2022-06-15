use std::{
    collections::{HashMap, BTreeMap},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    fs::File, error::Error,
};

use crate::{io::vcf::reader::vcfreader::VCFReader};
use crate::io::vcf::sampletag::SampleTag;

use log::{warn};
use rand::seq::SliceRandom;

/// Input sample panel definition file reader (`.panel` extension)
/// ### File characteristics:
/// - Tab-separated fields
/// - Columns: <Sample-id>  <Super-pop-id>    <pop-id>    <gender (optional)>
/// ### Fields:
/// - source_file: path to the source `.panel` file
/// - samples: HashMap with K: (String) <pop-id> | V: Vec<SampleTag>
#[derive(Debug, Clone)]
pub struct VCFPanelReader<'a> {
    pub source_file : &'a Path,
    pub samples: HashMap<String, Vec<SampleTag>>,
}



impl<'a> VCFPanelReader<'a> {
    /// Instantiate a Panel-definition from a `.panel` file
    /// # Arguments:
    /// - `panel_path`: path leading to the `.panel` file.
    /// # TODO: Convert this into ::from()
    pub fn new(panel_path: &Path) -> std::io::Result<VCFPanelReader> {
        let source = BufReader::new(File::open(panel_path)?);
        let samples = Self::parse(source)?;
        Ok(VCFPanelReader{samples, source_file: panel_path})
    }

    /// Copy and subset the source `.panel` files.
    /// The output file will only contain entries that are contained within the `samples` field.
    /// # Fields
    /// - `output_dir`: path leading to the output `.panel` file.
    pub fn copy_from_source(&self, output_dir: &Path) -> std::io::Result<()> {
        let source     = BufReader::new(File::open(self.source_file)?);
        let mut writer = BufWriter::new(File::create(output_dir)?);
        for line in source.lines(){
            let line = line?;
            let split_line = line.split('\t').collect::<Vec<&str>>();
            let (id, pop) = (split_line[0], split_line[1]);
            if self.samples[pop].contains(&SampleTag::new(id, None)) {
                writer.write(format!("{line}\n").as_bytes())?;
            }
        }
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
    pub fn fetch_contaminants(&self, contam_pop: &Vec<String>, contam_num_ind: &Vec<usize>) -> Result<Vec<Vec<SampleTag>>, Box<dyn Error>> {
        let mut samples_contam_tags: Vec<Vec<SampleTag>> = Vec::new();

        // ---- Loop along `contam_num_ind` and populate `samples_contam_tags` with n sampletags.
        for num in contam_num_ind {
            let mut contam_ind_tags : Vec<SampleTag> = Vec::new();

            // ---- Wrap if contam_pop.len() < contam_num_ind.len() OR if contam_pop.len() > contam_num_ind.len()
            let contam_pop = &contam_pop[*num % std::cmp::min(contam_num_ind.len(), contam_pop.len())];

            // ---- Randomly sample n sampletags
            for _ in 0..*num {
                let contam_sample_tag = self.random_sample(contam_pop, None)?.unwrap().clone();
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
    pub fn random_sample(&self, pop: &String, exclude: Option<&Vec<&SampleTag>>) -> Result<Option<&SampleTag>, Box<dyn Error>> {
        // Extract the vector of candidate SampleTags using the provided population-id.
        // Bailout if `pop` does not match anything.
        let candidates = match self.samples.get(pop) {
            Some(sample_vec) => sample_vec,
            None => return Err(format!("Could not fetch random sample using population tag: \"{pop}\"").into())
        };

        let candidate = match exclude {
            // If there are sample tags to exclude, pre-filter out these SampleTags from `self.samples` before random sampling.
            Some(tag_vec) => {
                candidates.iter()
                    .filter(|sample| ! tag_vec.contains(&sample))
                    .collect::<Vec<&SampleTag>>()
                    .choose(&mut rand::thread_rng()).map(|tag| *tag)
            },
            // If there are no sampleTag to exclude, directly proceed to random sampling.
            None => candidates.choose(&mut rand::thread_rng())
        };
        Ok(candidate)
    }

    /// Subset `self.samples` with entries matching the provided `subset_pops`.
    /// # Arguments:
    /// - `subset_pops`: slice of (super-)population-id.
    pub fn subset_panel(&mut self, subset_pops: &[String]) {
        self.samples.retain(|key,_| subset_pops.contains(key));
    }

    /// Convert the `self.samples` HashMap<String, Vec<SampleTag>> into a transposed BTreeMap<SampleTag, String>.
    pub fn into_transposed_btreemap(&self) -> BTreeMap<SampleTag, String> {
        let mut transposed_data = BTreeMap::new();
        for (key, values) in self.samples.iter(){
            for value in values.iter() {
                transposed_data.insert(value.clone(), key.clone());
            }
        }
        transposed_data
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
            output.entry(line[1].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx)); // Key == Super-pop
            output.entry(line[2].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx)); // Key == Pop
        }
        Ok(output)
    }

    /// Read a VCF file, search through its header for samples indices, and assign a field-index to each SampleTag found within `self.samples`
    /// # Arguments:
    /// - `vcf`: Path to the  targeted `.vcf` file. 
    pub fn assign_vcf_indexes(&mut self, vcf: &Path)  -> std::io::Result<()> {

        let reader = VCFReader::new(vcf, 0)?;
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
