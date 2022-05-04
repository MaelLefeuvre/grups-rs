use std::{
    collections::{HashMap, BTreeMap},
    io::{BufRead, BufReader},
    path::Path,
    fs::File,
};

use crate::io::vcf::reader::vcfreader::VCFReader;
use crate::io::vcf::sampletag::SampleTag;

use log::{warn};
use rand::seq::SliceRandom;

// TODO: Convert this into named tuple with deref.
#[derive(Debug, Clone)]
pub struct VCFPanelReader{
    pub samples: HashMap<String, Vec<SampleTag>>,
}



impl VCFPanelReader {
    // TODO: Convert this into ::from()
    pub fn new(panel_path: &Path) -> std::io::Result<VCFPanelReader> {
        let source = BufReader::new(File::open(panel_path)?);
        //let reader = VCFReader::new(vcf, 0)?;
        let samples = Self::parse_header(source)?;
        //drop(reader);
        Ok(VCFPanelReader{samples})
    }

    pub fn random_sample(&self, pop: &String) -> Option<&SampleTag> {
        self.samples[pop].choose(&mut rand::thread_rng())
    }

    pub fn subset_panel(&mut self, subset_pops: &Vec<String>) {
        self.samples.retain(|key,_| subset_pops.contains(key));
    }

    pub fn into_transposed_btreemap(&self) -> BTreeMap<SampleTag, String> {
        let mut transposed_data = BTreeMap::new();
        for (key, values) in self.samples.iter(){
            for value in values.into_iter() {
                transposed_data.insert(value.clone(), key.clone());
            }
        }
        transposed_data
    }

    pub fn sort_panel(&mut self){
        self.samples.iter_mut().for_each(|(_, tag_vec)| tag_vec.sort());
    }

    pub fn parse_header(source: BufReader<File>) -> std::io::Result<HashMap<String, Vec<SampleTag>>> {
        let mut output: HashMap<String, Vec<SampleTag>> = HashMap::new();

        for line in source.lines(){
            let line = line?;
            let line: Vec<&str> = line.split('\t').collect();
            let sample_idx = None;
            output.entry(line[1].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx));
            output.entry(line[2].into()).or_insert(Vec::new()).push(SampleTag::new(line[0], sample_idx));

        }
        Ok(output)
    }

    pub fn assign_vcf_indexes(&mut self, vcf: &Path)  -> std::io::Result<()> {

        let reader = VCFReader::new(vcf, 0)?;
        let header = &reader.samples();

        for (pop, sample_tags) in self.samples.iter_mut(){
            for sample_tag in sample_tags.iter_mut() {
                let sample_idx = match header.iter().position(|id| id == sample_tag.id()){
                    Some(idx) => idx - 9, // Correct for previous fields.
                    None => {warn!("Sample not found in input vcfs. Skipping:\n{:?}", pop); continue}
                };
                sample_tag.set_idx(sample_idx);
            }
        }
        Ok(())
    }
}
