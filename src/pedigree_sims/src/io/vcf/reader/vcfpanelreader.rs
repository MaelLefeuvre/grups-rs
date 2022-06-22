use std::{
    collections::{HashMap, BTreeMap},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    fs::File,
};

use crate::{io::vcf::reader::vcfreader::VCFReader};
use crate::io::vcf::sampletag::SampleTag;

use log::{warn};
use rand::seq::SliceRandom;

// TODO: Convert this into named tuple with deref.
#[derive(Debug, Clone)]
pub struct VCFPanelReader<'a> {
    pub source_file : &'a Path,
    pub samples: HashMap<String, Vec<SampleTag>>,
}



impl<'a> VCFPanelReader<'a> {
    // TODO: Convert this into ::from()
    pub fn new(panel_path: &Path) -> std::io::Result<VCFPanelReader> {
        let source = BufReader::new(File::open(panel_path)?);
        let samples = Self::parse_header(source)?;
        Ok(VCFPanelReader{samples, source_file: panel_path})
    }

    pub fn copy_from_source(&self, output_dir: &Path) -> std::io::Result<()> {

        println!("source: {:?}", self.source_file);
        println!("output: {:?}", output_dir);


        let source = BufReader::new(File::open(self.source_file)?);
        let mut writer = BufWriter::new(File::create(output_dir)?);
        for line in source.lines(){
            let line = line?;
            let split_line = line.split('\t').collect::<Vec<&str>>();
            let (id, pop) = (split_line[0], split_line[1]);
            if self.samples.get(pop).unwrap_or(&Vec::new()).contains(&SampleTag::new(id, None)) {
                writer.write(format!("{line}\n").as_bytes())?;
            }
        }
        Ok(())
    }

    pub fn fetch_contaminants(&self, contam_pop: &Vec<String>, contam_num_ind: &Vec<usize>) -> Vec<Vec<SampleTag>> {
        let mut samples_contam_tags: Vec<Vec<SampleTag>> = Vec::new();

        for num in contam_num_ind {
            let mut contam_ind_tags : Vec<SampleTag> = Vec::new();
            let contam_pop = &contam_pop[*num % std::cmp::min(contam_num_ind.len(), contam_pop.len())]; // Wrap if contam_pop.len() < contam_num_ind.len() OR if contam_pop.len() > contam_num_ind.len()

            for _ in 0..*num {
                let contam_sample_tag = self.random_sample(contam_pop, None).unwrap().clone();
                contam_ind_tags.push(contam_sample_tag);
            }
            samples_contam_tags.push(contam_ind_tags);
        }
        samples_contam_tags
    }

    pub fn random_sample(&self, pop: &String, exclude: Option<&Vec<&SampleTag>>) -> Option<&SampleTag> {
        let candidate = match exclude {
            Some(tag_vec) => {
                self.samples[pop].iter()
                    .filter(|sample| ! tag_vec.contains(&sample))
                    .collect::<Vec<&SampleTag>>()
                    .choose(&mut rand::thread_rng()).map(|tag| *tag)
            }
            None => self.samples[pop].choose(&mut rand::thread_rng())
        };
        candidate
    }

    pub fn subset_panel(&mut self, subset_pops: &[String]) {
        self.samples.retain(|key,_| subset_pops.contains(key));
    }

    pub fn into_transposed_btreemap(&self) -> BTreeMap<SampleTag, String> {
        let mut transposed_data = BTreeMap::new();
        for (key, values) in self.samples.iter(){
            for value in values.iter() {
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
