use std::{ops::Deref, sync::Arc, fmt::{self, Formatter, Display}};

use parking_lot::RwLock;

use super::Individual;

/// Helper struct to Deref and Display the parents of an individual.
/// 
#[derive(Debug, Clone)]
pub struct Parents([Arc<RwLock<Individual>>; 2]);

impl Parents{
    /// Instantiate a new set of Parents.
    /// # Arguments
    /// - `parents`: size-two array of Arc<RwLock<Individual>. Each Individual being a parent.  
    pub fn new(parents: [Arc<RwLock<Individual>>; 2]) -> Parents {
        Parents(parents)
    }
}

impl PartialEq for Parents {
    fn eq(&self, other: &Self) -> bool {
        *self.0[0].read() == *other.0[0].read() && 
        *self.0[1].read() == *other.0[1].read()
    }
}

impl Deref for Parents {
    type Target = [Arc<RwLock<Individual>>];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Display for Parents {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{} <-> {}", 
            &self[0].read().label,
            &self[1].read().label
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
        let father = Arc::new(RwLock::new(common::mock_founder(parents_labels[0])));
        let mother = Arc::new(RwLock::new(common::mock_founder(parents_labels[1])));
        let parents = Parents::new([father, mother]);
        let display = format!("{parents}");

        assert!(display.contains("father"));
        assert!(display.contains("mother"));
    }
}