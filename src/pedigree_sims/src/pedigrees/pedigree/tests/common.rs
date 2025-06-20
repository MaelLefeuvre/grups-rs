use crate::pedigrees::pedigree::Individual;
use crate::pedigrees::pedigree::PedComparison;

use std::sync::Arc;
use parking_lot::RwLock;

/// Mock a founder individual with no parents.
/// # Arguments:
///  - `label`: raw string slice defining the individual's name (e.g. "father", "mother", etc.)
pub fn mock_founder(label: &str) -> Individual {
    Individual::new(label, None, None)
}

/// Mock an offspring individual, with two parents.
/// # Arguments:
///  - `label` : raw string slice defining the individual's name (e.g. "child")
/// - `parents`: sice two array, defining the parents names (e.g. `["mother", "father"]``)
pub fn mock_offspring(label: &str, parents_labels: Option<[&str; 2]>) -> Individual {
    let parents_labels = parents_labels.unwrap_or(["father", "mother"]);

    let father = Arc::new(RwLock::new(mock_founder(parents_labels[0])));
    let mother = Arc::new(RwLock::new(mock_founder(parents_labels[1])));
    let parents = Some([&father, &mother].into());
    Individual::new(label, parents, None)
}

/// Mock a pedigree comparison, using two founder individuals (=> this comparison is always labeled "unrelated")
pub fn mock_pedcomparison() -> PedComparison {
    let ind1 = Arc::new(RwLock::new(mock_founder("ind1")));
    let ind2 = Arc::new(RwLock::new(mock_founder("ind2")));
    PedComparison::new("unrelated", [&ind1, &ind2], false)
}