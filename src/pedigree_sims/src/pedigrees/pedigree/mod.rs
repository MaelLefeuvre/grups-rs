mod individual;
use individual::Individual;

mod pedigree;
pub use pedigree::Pedigree;

pub mod io;

mod contaminant;
pub use contaminant::Contaminant;

mod comparisons;
use comparisons::PedComparisons;
use comparisons::PedComparison;


pub mod pedigree_params;
use pedigree_params::PedigreeParams;
