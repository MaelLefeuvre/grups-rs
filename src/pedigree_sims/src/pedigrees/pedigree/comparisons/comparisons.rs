use std::ops::{Deref, DerefMut};

use super::PedComparison;

#[derive(Debug, Clone)]
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

impl std::fmt::Display for PedComparisons {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0.iter().fold(Ok(()), |result, comparison| {
            result.and_then(|_| writeln!(f, "{}", comparison))
        })
    }
}

impl PedComparisons {
    pub fn new() -> Self {
        Self(Vec::new())
    }
}