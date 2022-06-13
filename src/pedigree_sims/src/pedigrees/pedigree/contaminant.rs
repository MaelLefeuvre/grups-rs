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

impl std::fmt::Display for Contaminant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.iter().fold(Ok(()), |_, contaminants| {
            write!(f, "- [")?;
            contaminants.iter().fold(Ok(()), |_, contaminant| {
                write!(f, " {contaminant} ")
            })?;
            writeln!(f, "]")
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::genotype_reader::MockGenotypeReader;

    fn dummy_sample_tags(nums: [usize; 2]) -> [Vec<SampleTag>; 2] {
        let mut sample_tags = [vec![], vec![]];
        let mut id = 0;
        for (i, num) in nums.into_iter().enumerate() {
            for _ in 0..num {
                let new_tag = SampleTag::new(id.to_string().as_str(), None);
                sample_tags[i].push(new_tag);
                id+=1;
            }
        }
        sample_tags
    }


    #[test]
    fn flat_list_getter() {
        for i in 0..10 {
            for j in 0..10 {
                let tag_numbers = [i,j];
                let contaminant = Contaminant::new(dummy_sample_tags(tag_numbers));
                let flat_list = contaminant.as_flat_list();
                let expected_len: usize = tag_numbers.iter().sum();
                assert_eq!(flat_list.len(), expected_len);
            }
        }
    }

    #[test]
    fn local_cont_af(){
        let tag_numbers = [2,2];
        let contaminant = Contaminant::new(dummy_sample_tags(tag_numbers));
        let mut mock_reader = MockGenotypeReader::default();

        let mut dummy_alleles = vec![
            Some([0,0]), // Alleles of ind 2 - sampletag 2 | -> exp_af = 0.25
            Some([0,1]), // Alleles of ind 2 - sampletag 1 | 
            Some([1,0]), // Alleles of ind 1 - sampletag 2 | -> exp_af = 0.75
            Some([1,1])  // Alleles of ind 1 - sampletag 1 | 
        ]; 

        mock_reader.expect_get_alleles() // Expected output -> [0.75, 0.25]
            .times(dummy_alleles.len())
            .returning(move |_| {
                dummy_alleles.pop().unwrap() // pop() -> remember we're iterating in reverse.
            });

        let want = [0.75, 0.25];
        let got = contaminant.compute_local_cont_af(&mock_reader).unwrap();

        assert_eq!(want, got)
    }
}