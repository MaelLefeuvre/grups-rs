use crate::snpcoord::SNPCoord;
use crate::genome::Genome;
use std::collections::HashMap;
use itertools::Itertools;
use std::fmt;

use super::JackknifeBlock;

#[derive(Debug)]
/// A simple struct representing a HashMap of Jackknife Blocks for a given
/// genome. Implemented within struct Comparison. See: `pwd_from_stdin::pileup::Comparison`
///   - HashMap is indexed according to the chromosome name.
/// 
/// # Traits:
///   - `Display` : Pretty print. for file writing and debug. Recursively calls `Display` for each
///                 `JackknifeBlock` found within the `blocks` field.
/// # TODO:
///   - add a header for the `Display trait`
pub struct JackknifeBlocks {
    blocks : HashMap<u8, Vec<JackknifeBlock>>
}

impl JackknifeBlocks {
    pub fn new(genome: &Genome, blocksize: u32) -> JackknifeBlocks {
        let mut jackknives = HashMap::new();
        for chr in genome.values() {
            let mut blocks = Vec::new();
            let mut add_block = |chr, start, end| {
                blocks.push(JackknifeBlock::new(chr, start, end))
            };
            let chr_blocks: Vec<u32> = (1..chr.length+1).step_by(blocksize as usize).collect();
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

    /// Search for a given block, using an SNPCoord struct.
    /// Return the block which contains the SNPCoord position. 
    pub fn find_block(&mut self, coordinate: &SNPCoord) -> Option<&mut JackknifeBlock> {
        self.blocks.get_mut(&coordinate.chromosome)
            .unwrap()
            .iter_mut()
            .find(|block| block.range.contains(&coordinate.position))
    }
}

// Good Stuff: https://github.com/apolitical/impl-display-for-vec
impl fmt::Display for JackknifeBlocks {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.blocks.keys().sorted().fold(Ok(()), |_result, chr  | {
            self.blocks[chr].iter().fold(Ok(()), |result, block| {
                result.and_then(|_| writeln!(f, "{}", block))
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
       let genome = Genome::from(&vec![Chromosome::new(0, 1, 249250621)]);
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[&1].name].len(),249251);
    }
    
    #[test]
    fn jackknife_blocks_check_len_equal () {
       let genome = Genome::from(&vec![Chromosome::new(0, 1, 2000001)]);
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[&1].name].len(), 2000);
    }

    #[test]
    fn display(){
        let genome = Genome::from(&vec![ 
            Chromosome::new(0, 1, 10_001),
            Chromosome::new(1, 2, 10_001), 
        ]);
        
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
                ).unwrap();
            }
        }
        assert_eq!(expected_output, format!("{blocks}"));
    }

    #[test]
    fn search_block() {
        let genome = Genome::default();
        let blocksize = 1000;
        let mut blocks = JackknifeBlocks::new(&genome, blocksize);

        let coord = SNPCoord{chromosome: 1, position: 651_000, reference: None, alternate: None};

        let block = blocks.find_block(&coord).unwrap();
        assert_eq!(block.chromosome, coord.chromosome);
        assert!(block.range.start <= coord.position);
        assert!(block.range.end > coord.position);
    }
}

