pub mod coordinate;

pub mod nucleotide;
pub use nucleotide::Nucleotide;
pub use nucleotide::Phred;

pub mod snp;
pub use snp::SNPCoord;

pub mod jackknife;

mod genetic_map;
pub use genetic_map::GeneticMap;

mod genome;
pub use crate::genome::Genome; 

mod chromosome;
pub use chromosome::Chromosome;


