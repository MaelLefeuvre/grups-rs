pub mod coordinate;

pub mod snp;
pub use snp::SNPCoord;

pub mod jackknife;

mod genetic_map;
pub use genetic_map::GeneticMap;

mod genome;
pub use crate::genome::Genome; 

mod chromosome;
pub use crate::chromosome::Chromosome;


