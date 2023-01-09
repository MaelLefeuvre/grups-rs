mod individual;
pub use individual::Individual;

mod comparison;
#[allow(clippy::module_inception)]
mod comparisons;
pub use comparisons::Comparisons;

pub mod pwd;
pub use pwd::Pwd;

#[cfg(test)]
mod tests;

const UNDEFINED_LABEL_PREFIX: &str = "Ind";
const DISPL_SEP             : &str  = " - ";
const PAIRS_FORMAT_LEN      : usize = 20;
const COUNT_FORMAT_LEN      : usize = 11;
const AVERG_FORMAT_LEN      : usize = 11;
const FLOAT_FORMAT_PRECISION: usize = 6;