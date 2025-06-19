use crate::genome::Genome;
use itertools::Itertools;
use std::fmt::{self, Display, Formatter};

use super::JackknifeBlock;
use crate::coordinate::{ChrIdx, Coordinate};

use ahash::AHashMap;

#[derive(Debug)]
pub struct JackknifeEstimates {
    pub estimate: f64,
    pub variance: f64,
}

#[derive(Debug)]
/// A simple struct representing a `HashMap` of `JackknifeBlock` for a given
/// genome. Implemented within struct Comparison. See: `pwd_from_stdin::pileup::Comparison`
///   - `HashMap` is indexed according to the chromosome name.
/// 
/// # Traits:
///   - `Display` : Pretty print. for file writing and debug. Recursively calls `Display` for each
///     `JackknifeBlock` found within the `blocks` field.
/// # TODO:
///   - add a header for the `Display trait`
pub struct JackknifeBlocks {
    blocks : AHashMap<ChrIdx, Vec<JackknifeBlock>>
}

impl JackknifeBlocks {
    #[must_use]
    pub fn new(genome: &Genome, blocksize: u32) -> JackknifeBlocks {
        let mut jackknives = AHashMap::with_capacity(genome.len());
        for chr in genome.values() {
            let mut blocks = Vec::with_capacity((chr.length / blocksize) as usize);
            let mut add_block = |chr, start, end| {
                blocks.push(JackknifeBlock::new(chr, start, end));
            };
            let chr_blocks: Vec<u32> = (1..=chr.length).step_by(blocksize as usize).collect();
            for i in 1..chr_blocks.len(){
                add_block(chr.name, chr_blocks[i-1], chr_blocks[i]);
            }
            if (chr.length) % (blocksize) != 1 {
                add_block(chr.name, chr_blocks[chr_blocks.len()-1], chr.length);
            }
            jackknives.insert(chr.name, blocks);
        }
        JackknifeBlocks{blocks: jackknives}
    }

    /// Search for a given block, using an `SNPCoord` struct.
    /// Return the block containing the `SNPCoord` position. 
    pub fn find_block(&mut self, coordinate: &Coordinate) -> Option<&mut JackknifeBlock> {
        self.blocks.get_mut(&coordinate.chromosome)?
            .iter_mut()
            .find(|block| block.range.contains(&coordinate.position))
    }

    /// Compute jackknifed avg. PWD estimate for each chromosome block
    #[must_use]
    pub fn compute_unequal_delete_m_pseudo_values(&self, sum_pwd: f64, sum_overlap: u32) -> JackknifeEstimates {
        let mut theta_jk: f64 = 0.0;
        for chromosome_blocks in self.blocks.values() {
            for block in chromosome_blocks.iter() {
                let pseudo_value = block.compute_unequal_delete_m_pseudo_value(sum_pwd, sum_overlap);
                if pseudo_value.hj.is_finite() {
                    theta_jk += pseudo_value.weigthed_pseudovalue();
                };
            }
        }

        // Compute Jackknife variance estimate 
        let mut var_jk  : f64 = 0.0;
        for chromosome_blocks in self.blocks.values() {
            for block in chromosome_blocks.iter() {
                let pseudo_value = block.compute_unequal_delete_m_pseudo_value(sum_pwd, sum_overlap);
                if pseudo_value.hj.is_finite() {
                    var_jk += f64::powf(pseudo_value.weigthed_pseudovalue() - theta_jk, 2.0) / (pseudo_value.hj - 1.0);
                }
            }
        }

        let g = self.blocks.values().map(Vec::len).sum::<usize>() as u32; // Casting to u32 because usize -> f64 conversion
        let var_jk = var_jk / g as f64;                                   // can generate precision loss on x64 architectures.

        JackknifeEstimates{estimate: theta_jk, variance: var_jk}
    }

    
}

// Good Stuff: https://github.com/apolitical/impl-display-for-vec
impl Display for JackknifeBlocks {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        self.blocks.keys().sorted().try_fold((), |_result, chr  | {
            self.blocks[chr].iter().try_fold((), |_, block| {
                writeln!(f, "{block}")
            })
        })
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    use std::fmt::Write;
    use super::super::{CHROM_FORMAT_LEN, COUNT_FORMAT_LEN, RANGE_FORMAT_LEN};
    use crate::chromosome::Chromosome;

    #[test]
    fn jackknife_blocks_check_len_unequal() {
       let genome = Genome::from(&[Chromosome::new(1, 249250621)]);
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[&ChrIdx::from(1)].name].len(),249251);
    }
    
    #[test]
    fn jackknife_blocks_check_len_equal () {
       let genome = Genome::from(&[Chromosome::new(1, 2000001)]);
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[&ChrIdx::from(1)].name].len(), 2000);
    }

    #[test]
    fn display(){
        let genome = Genome::from(&[
            Chromosome::new(1, 10_001),
            Chromosome::new(2, 10_001)]);
        
        let blocksize = 1000;
        let blocks = JackknifeBlocks::new(&genome, blocksize);

        // Expected output
        let mut expected_output = String::new();
        for chr in genome.values() {
            let expected_ranges: Vec<u32> = (1..chr.length+1).step_by(blocksize as usize).collect();
            for step in expected_ranges.windows(2)  {
                let (start, end) = (step[0], step[1]);

                writeln!(expected_output,
                    "{: <CHROM_FORMAT_LEN$} - \
                    {: <RANGE_FORMAT_LEN$} - \
                    {: <RANGE_FORMAT_LEN$} - \
                    {: <COUNT_FORMAT_LEN$} - \
                    {: <COUNT_FORMAT_LEN$}",
                    chr.name, start, end, 0, 0
                ).expect("Jackknife Blocks 'display()' test error while writing line");
            }
        }
        assert_eq!(expected_output, format!("{blocks}"));
    }

    #[test]
    fn search_block() {
        let genome = Genome::default();
        let blocksize = 1000;
        let mut blocks = JackknifeBlocks::new(&genome, blocksize);

        let coord = Coordinate::new(1, 651_000);

        let block = blocks
            .find_block(&coord)
            .expect("Jackknife Blocks search_block() test error. Failed to find block");
        assert_eq!(block.chromosome, coord.chromosome);
        assert!(block.range.start <= coord.position);
        assert!(block.range.end > coord.position);
    }
}

