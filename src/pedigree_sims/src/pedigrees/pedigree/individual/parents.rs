use std::{
    ops::Deref,
    cell::RefCell,
    rc::Rc,
};

use super::Individual;

/// Helper struct to Deref and Display the parents of an individual.
/// 
#[derive(Debug, Clone, PartialEq)]
pub struct Parents([Rc<RefCell<Individual>>; 2]);

impl Parents{
    /// Instantiate a new set of Parents.
    /// # Arguments
    /// - `parents`: size-two array of Rc<RefCell<Individual>. Each Individual being a parent.  
    pub fn new(parents: [Rc<RefCell<Individual>>; 2]) -> Parents {
        Parents(parents)
    }
}

impl Deref for Parents {
    type Target = [Rc<RefCell<Individual>>; 2];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::fmt::Display for Parents {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} <-> {}", 
            RefCell::borrow(&self[0]).label,
            RefCell::borrow(&self[1]).label
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
        let father = Rc::new(RefCell::new(common::mock_founder(parents_labels[0])));
        let mother = Rc::new(RefCell::new(common::mock_founder(parents_labels[1])));
        let parents = Parents::new([father, mother]);
        let display = format!("{}", parents);

        assert!(display.contains("father"));
        assert!(display.contains("mother"));
    }
}