use std::ops::Deref;

use located_error::prelude::*;

mod error;
use error::InfoFieldError;

/// Struct representing the VCF `INFO` field.
#[derive(Debug, Default)]
pub struct InfoField (Option<Vec<String>>);

impl Deref for InfoField {
    type Target = Vec<String>;
    fn deref(&self) -> &Self::Target {
        self.0.as_ref().expect("Attempting to access empty InfoField!")
    }
}

impl InfoField {
    /// Instantiate a new InfoField from a provided string slice.
    /// # Arguments:
    /// - `field`: raw, unparsed `INFO` vcf column string slice.
    pub fn new(field: &str) -> Self {
        Self(Some(field.split(';').map(|s| s.to_string()).collect::<Vec<String>>()))
    }

    /// Clear-out the structs content.
    pub fn clear(&mut self) {
        self.0 = None
    }

    /// return `true` if InfoField contains the "MULTI_ALLELIC" tag. (i.e. is multiallelic.) 
    pub fn is_multiallelic(&self) -> bool {
        self.iter().any(|field| field == "MULTI_ALLELIC")
    }

    /// return `true` if InfoField contains the "VT=SNP" tag. (i.e. is an SNP)
    pub fn is_snp(&self) -> Result<bool> {
        let vtype = self.iter()
            .find(|&field| field.starts_with("VT="))
            .ok_or(InfoFieldError::MissingVT)?
            .split('=')
            .collect::<Vec<&str>>()[1];
        Ok(vtype == "SNP")
    }

    /// Return the annotated population allele frequency for a given population.
    /// # Arguments:
    /// - `pop`: super-population-id. 
    /// 
    /// # Behavior: 
    /// for a given `pop` raw string slice, this method will search for any field matching "{pop}_AF="
    /// and return its value.
    pub fn get_pop_allele_frequency(&self, pop: &str) -> Result<f64> {
        let pop_af_regex = format!("{}_AF", pop);
        self.iter()
            .find(|&field| field.starts_with(&pop_af_regex))
            .and_then(|x| x.split('=').last())
            .ok_or(InfoFieldError::MissingAF(pop_af_regex))?
            .parse::<f64>()
            .map_err(InfoFieldError::ParseAlleleFrequency)
            .loc("While attempting to retrieve the allele frequency of population '{pop}'")
    }
}