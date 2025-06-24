
use crate::pedigrees::pedigree::Pedigree;

use super::super::{
    Individual, IndividualId,
    Relationship, RelationshipId,
    PedComparison,
};

use std::hash::{DefaultHasher, Hash, Hasher};

/// Mock a founder individual with no parents.
/// # Arguments:
///  - `label`: raw string slice defining the individual's name (e.g. "father", "mother", etc.)
pub fn mock_founder_pedigree(label: &str) -> Pedigree {
    //let mut hasher = DefaultHasher::new();
    //label.hash(&mut hasher);
    //let hashed = hasher.finish().into();
    //let mock_id = IndividualId::from_u64(hashed);
    //Individual::new(mock_id, label, None)

    let mut pedigree = Pedigree::new();
    pedigree.add_individual(label, None, None);
    //pedigree.get_ind(iid).unwrap().clone()
    pedigree
}

/// Mock an offspring individual, with two parents.
/// # Arguments:
///  - `label` : raw string slice defining the individual's name (e.g. "child")
/// - `parents`: sice two array, defining the parents names (e.g. `["mother", "father"]``)
pub fn mock_offspring_pedigree(label: &str, parents_labels: Option<[&str; 2]>) -> Pedigree {
    let parents_labels = parents_labels.unwrap_or(["father", "mother"]);

    //let mut child = mock_founder(label);
    //let relationships = [
    //    RelationshipId::from_u64(0),
    //    RelationshipId::from_u64(1),
    //];
    //child.add_parents(relationships);
    //child

    let mut pedigree = Pedigree::new();
    let _fid = pedigree.add_individual(parents_labels[0], None, None);
    let _mid = pedigree.add_individual(parents_labels[1], None, None);
    let _cid = pedigree.add_individual(label, Some(parents_labels), None);
    //pedigree.set_relationship(cid, [fid, mid]);
    pedigree
}

/// Mock a pedigree comparison, using two founder individuals (=> this comparison is always labeled "unrelated")
pub fn mock_pedcomparison() -> Pedigree {
    //let ind1 = Arc::new(RwLock::new(mock_founder("ind1")));
    //let ind2 = Arc::new(RwLock::new(mock_founder("ind2")));
    //let ind1 = mock_founder("ind1").id;
    //let ind2 = mock_founder("ind2").id;
    //PedComparison::new("unrelated", [ind1, ind2], false)

    let mut pedigree = Pedigree::new();
    let _fid = pedigree.add_individual("ind1", None, None);
    let _mid = pedigree.add_individual("ind2", None, None);
    pedigree.add_comparison("unrelated", ["ind1", "ind2"]).unwrap();
    //let comp = pedigree.comparisons.get(comp_id).unwrap().clone()
    pedigree
}