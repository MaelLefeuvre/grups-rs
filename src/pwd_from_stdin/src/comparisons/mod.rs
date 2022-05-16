mod individual;
pub use individual::Individual;

mod comparison;

mod comparisons;
pub use comparisons::Comparisons;

mod tests;

const UNDEFINED_LABEL_PREFIX: &str = "Ind";
const DISPL_SEP             : &str  = " - ";
const PAIRS_FORMAT_LEN      : usize = 20;
const COUNT_FORMAT_LEN      : usize = 7;
const AVERG_FORMAT_LEN      : usize = 8;
const FLOAT_FORMAT_PRECISION: usize = 5;