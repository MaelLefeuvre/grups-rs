
use crate::pedigrees::pedigree::Pedigree;

/// Mock a founder individual with no parents.
/// # Arguments:
///  - `label`: raw string slice defining the individual's name (e.g. "father", "mother", etc.)
pub fn mock_founder_pedigree(label: &str) -> Pedigree {
    let mut pedigree = Pedigree::new();
    pedigree.add_individual(label, None, None);
    pedigree
}

/// Mock an offspring individual, with two parents.
/// # Arguments:
///  - `label` : raw string slice defining the individual's name (e.g. "child")
/// - `parents`: sice two array, defining the parents names (e.g. `["mother", "father"]``)
pub fn mock_offspring_pedigree(label: &str, parents_labels: Option<[&str; 2]>) -> Pedigree {
    let parents_labels = parents_labels.unwrap_or(["father", "mother"]);

    let mut pedigree = Pedigree::new();
    let _fid = pedigree.add_individual(parents_labels[0], None, None);
    let _mid = pedigree.add_individual(parents_labels[1], None, None);
    let _cid = pedigree.add_individual(label, Some(parents_labels), None);
    pedigree
}

/// Mock a pedigree comparison, using two founder individuals (=> this comparison is always labeled "unrelated")
pub fn mock_pedcomparison() -> Pedigree {
    let mut pedigree = Pedigree::new();
    let _fid = pedigree.add_individual("ind1", None, None);
    let _mid = pedigree.add_individual("ind2", None, None);
    pedigree.add_comparison("unrelated", ["ind1", "ind2"]).expect("Individual should be includable");
    pedigree
}