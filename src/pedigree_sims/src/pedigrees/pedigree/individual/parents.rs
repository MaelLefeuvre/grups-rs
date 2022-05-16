use std::{
    ops::Deref,
    cell::RefCell,
    rc::Rc,
};

use super::Individual;



#[derive(Debug, Clone, PartialEq)]
pub struct Parents([Rc<RefCell<Individual>>; 2]);

impl Parents{
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