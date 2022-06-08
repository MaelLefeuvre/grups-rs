
const CHROM_FORMAT_LEN: usize = 4;
const RANGE_FORMAT_LEN: usize = 15;
const COUNT_FORMAT_LEN: usize = 5;

mod jackknife_block;
pub use jackknife_block::JackknifeBlock;

mod jackknife_blocks;
pub use jackknife_blocks::JackknifeBlocks;
pub use jackknife_blocks::JackknifeEstimates;