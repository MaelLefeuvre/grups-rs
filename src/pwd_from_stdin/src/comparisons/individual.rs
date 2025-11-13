
use crate::pileup::Pileup;


use super::{UNDEFINED_LABEL_PREFIX, error::ComparisonError};

/// Represents a requested individual Within the pileup.
///  - `name`      : Name of the individual. Either given through user-input, or constructed as `Ind{index}` by default
///    when no name has been provided.
/// 
///  - `index`     : 0 based index of the individual within the pileup. Note that this index is technically offset by 3,
///    since the first three columns of a pileup respectively define 'chr', 'pos', 'ref'.
/// 
///  - `min_depth` : minimum sequencing depth that is allowed before making a comparison. User-defined, or defaults to 1.
/// 
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Individual {
    pub name     : String,
    pub index    : usize,
    pub min_depth: u16,
}

impl Individual {
    /// Instantiate a new Pileup Individual, with a provided name, minimum requested depth and index.
    #[must_use]
    pub fn new (name: Option<&String>, index: usize, min_depth: u16) -> Individual {
        // Parse name and give default if non-existent
        let name = match name {                            
            Some(name) => name.to_string(),
            None => format!("{UNDEFINED_LABEL_PREFIX}{index}"),
        };
        Individual{name, index, min_depth}
    }

    /// Ensure the current pileup depth is greater than the requested minimum depth at this line.
    pub fn satisfiable_depth(&self, pileups: &[Pileup]) -> Result<bool, ComparisonError> {
        pileups.get(self.index)
            .map(|pileup| pileup.depth >= self.min_depth)
            .ok_or(ComparisonError::InvalidPileupIndex(self.index, self.name.clone()))
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use genome::snp::Allele;

    use super::*;
    #[test]
    fn new_undefined_name() {
        let index = 0;
        let ind = Individual::new(None, 0, 0);
        assert_eq!(ind.name, format!("{UNDEFINED_LABEL_PREFIX}{index}"));
    }

    #[test]
    fn new_defined_name() {
        let test_name = "Test_name_label-01".to_string();
        let ind = Individual::new(Some(&test_name) , 0, 0);
        assert_eq!(ind.name, test_name);
    }

    #[test]
    fn satisfiable_depth() -> Result<(), Box<dyn Error>>{
        #![allow(clippy::cast_possible_truncation)]
        let (bases, scores) = ("AATA", "JEEJ");
        let pileup = [Pileup::new(Allele::N, bases.len() as u16, bases, scores, true)?];

        // ind.min_depth < pileup.depth   -> true
        let ind = Individual::new(None, 0, 1);
        assert!(ind.satisfiable_depth(&pileup)?);

        // ind.min_depth == pileup.depth -> true
        let ind = Individual::new(None, 0, 4);
        assert!(ind.satisfiable_depth(&pileup)?);

        // ind.min_depth > pileup.depth  -> false
        let ind = Individual::new(None, 0, 8);
        assert!(!ind.satisfiable_depth(&pileup)?);
        Ok(())
    }
}


