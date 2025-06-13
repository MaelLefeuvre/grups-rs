use std::{
    ops::Deref, sync::{Arc, Mutex}
};

use super::Individual;

/// Helper struct to Deref and Display the parents of an individual.
/// 
#[derive(Debug, Clone)]
pub struct Parents([Arc<Mutex<Individual>>; 2]);

impl Parents{
    /// Instantiate a new set of Parents.
    /// # Arguments
    /// - `parents`: size-two array of Arc<Mutex<Individual>. Each Individual being a parent.  
    pub fn new(parents: [Arc<Mutex<Individual>>; 2]) -> Parents {
        Parents(parents)
    }
}

impl PartialEq for Parents {
    fn eq(&self, other: &Self) -> bool {
        *self.0[0].lock().unwrap() == *other.0[0].lock().unwrap() && 
        *self.0[1].lock().unwrap() == *other.0[1].lock().unwrap()
    }
}

impl Deref for Parents {
    type Target = [Arc<Mutex<Individual>>];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::fmt::Display for Parents {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} <-> {}", 
            &self[0].lock().unwrap().label,
            &self[1].lock().unwrap().label
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pedigrees::pedigree::tests::common;
    #[test]
    fn display(){
        let parents_labels = ["father", "mother"];
        let father = Arc::new(Mutex::new(common::mock_founder(parents_labels[0])));
        let mother = Arc::new(Mutex::new(common::mock_founder(parents_labels[1])));
        let parents = Parents::new([father, mother]);
        let display = format!("{parents}");

        assert!(display.contains("father"));
        assert!(display.contains("mother"));
    }
}