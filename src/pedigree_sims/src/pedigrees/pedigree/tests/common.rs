use crate::pedigrees::pedigree::Individual;
use crate::pedigrees::pedigree::PedComparison;

use std::{
    rc::Rc,
    cell::RefCell,
};

pub fn mock_founder(label: &str) -> Individual {
    Individual::new(label, None)
}

pub fn mock_offspring(label: &str, parents_labels: Option<[&str; 2]>) -> Individual {
    let parents_labels = parents_labels.unwrap_or(["father", "mother"]);

    let father = Rc::new(RefCell::new(mock_founder(parents_labels[0])));
    let mother = Rc::new(RefCell::new(mock_founder(parents_labels[1])));
    let parents = Some([&father, &mother]);
    Individual::new(label, parents)
}

pub fn mock_pedcomparison() -> PedComparison {
    let ind1 = Rc::new(RefCell::new(mock_founder("ind1")));
    let ind2 = Rc::new(RefCell::new(mock_founder("ind2")));
    PedComparison::new("unrelated", [&ind1, &ind2], false)
}