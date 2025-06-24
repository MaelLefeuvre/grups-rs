use std::{fmt::{self, Display, Formatter}, ops::{Deref, DerefMut}};

mod error;
pub use error::ComparisonError;

mod comparison;
pub use comparison::PedComparison;

// ------------------------------------------------------------------ //
impl Display for PedComparisons {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        self.0.iter().try_fold((), |(), comparison| writeln!(f, "{comparison:?}"))
    }
}

/// Vector of all the tracked pedigree comparisons for a given Pedigree
#[derive(Debug, Clone, Default)]
pub struct PedComparisons(Vec<PedComparison>);

impl Deref for PedComparisons {
    type Target = Vec<PedComparison>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for PedComparisons {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl PedComparisons {
    /// Instantiate a new set of pedigree comparisons.
    pub fn new() -> Self {
        Self::default()
    }
}