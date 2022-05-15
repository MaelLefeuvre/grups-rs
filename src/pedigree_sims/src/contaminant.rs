use crate::io::{vcf::SampleTag, genotype_reader::GenotypeReader};

use std::{ops::{Deref, DerefMut}, error::Error};

#[derive(Debug)]
pub struct Contaminant([Vec<SampleTag>; 2]);

impl Deref for Contaminant {
    type Target = [Vec<SampleTag>; 2];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Contaminant {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Contaminant {
    pub fn new(sample_tags: [Vec<SampleTag>; 2]) -> Self {
        Contaminant(sample_tags)
    }

    pub fn as_flat_list(&self) -> Vec<&SampleTag> {
        let mut flat_list = Vec::with_capacity(self[0].len() + self[1].len());
        for vec in self.iter() {
            flat_list.extend(vec);
        }
        flat_list
    }

    pub fn compute_local_cont_af(&self, reader: &dyn GenotypeReader) -> Result<[f64; 2], Box<dyn Error>> {
        let mut output = Vec::with_capacity(2);

        for contaminant in self.0.iter() {
            let mut ref_allele_count = 0;
            let mut alt_allele_count = 0;
            for tag in contaminant.iter() {
                reader.get_alleles(tag).unwrap().into_iter().for_each(|allele| {
                    match allele {
                        0          => ref_allele_count += 1,
                        1          => alt_allele_count += 1,
                        other => panic!("Contaminating individual is multiallelic: {other}")
                    }
                })
            }
            let af = (alt_allele_count as f64) /(alt_allele_count as f64 + ref_allele_count as f64);
            output.push(af);
        }
        Ok(output.try_into().unwrap_or_else(|v: Vec<f64>| panic!("Expected a Vec of length {} but it was {}", 2, v.len())))
    }
}