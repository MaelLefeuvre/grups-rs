use std::{collections::BTreeMap, fmt::{self, Display, Formatter}, ops::{Deref, DerefMut}};

mod error;
pub use error::ComparisonError;

use grups_io::read::SampleTag;
use pwd_from_stdin::comparisons::Comparison;
use slotmap::SlotMap;

mod comparison;
pub use comparison::{PedComparison, ComparisonId};

use crate::pedigrees::pedigree::individual::IndividualId;

// ------------------------------------------------------------------ //
// ---- PedComparisons
#[derive(Debug, Clone, Default)]
pub struct PedComparisons{
    pub inner: SlotMap<ComparisonId, PedComparison>,
    pub labels: BTreeMap<String, ComparisonId>,
}

//impl Deref for PedComparisons {
//    type Target = SlotMap<ComparisonId, PedComparison>;
//    fn deref(&self) -> &Self::Target {
//        &self.ord
//    }
//}

//impl DerefMut for PedComparisons {
//    fn deref_mut(&mut self) -> &mut Self::Target {
//        &mut self.inner
//    }
//}

impl Display for PedComparisons {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        self.ordered().try_fold((), |(), comparison| writeln!(f, "{comparison:?}"))
    }
}

impl PedComparisons {
    /// Instantiate a new set of pedigree comparisons.
    pub fn new() -> Self {
        Self::default()
    }

    pub fn insert_with_key(&mut self, label: &str, pair_ids: [IndividualId; 2]) -> ComparisonId {
        if let Some(id) = self.labels.get(label) {
            *id
        } else {
            let comp_id = self.inner.insert_with_key(|id| PedComparison::new(id, label, pair_ids));
            *self.labels.entry(label.to_string()).or_insert(comp_id)
        }
    }

    pub fn get(&self, id: ComparisonId) -> Option<&PedComparison> {
        self.inner.get(id)
    }

    pub fn ordered(&self) -> impl Iterator<Item = &PedComparison >  {
        self.labels.values().copied().map(|id| self.inner.get(id).unwrap())
    }
}