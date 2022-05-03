use crate::snpcoord::SNPCoord;
use crate::genome::{Genome};
use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::ops::Range;
use itertools::Itertools;
use std::fmt;

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
    pub fn find_block(&mut self, coordinate: &SNPCoord) -> &mut JackknifeBlock {
        self.blocks.get_mut(&coordinate.chromosome)
            .unwrap()
            .iter_mut()
            .find(|block| block.range.contains(&coordinate.position)).unwrap()
    }

    /// Deprecated (JackknifeBlocks now implements the Display trait) 
    /// print each block + a header. Mainly for debug. 
    pub fn _print(&self) {
        println!("{: <4} - {: <15} - {: <15} - {: <5} - {: <5}", "chr", "start", "end", "overlap", "pwd");
        for chr in self.blocks.keys().sorted() {
            self.blocks[chr].iter().for_each(|block| println!("{}", block));
        }
    }
}

// Good Stuff: https://github.com/apolitical/impl-display-for-vec
impl fmt::Display for JackknifeBlocks {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //println!("{: <4} - {: <15} - {: <15} - {: <5} - {: <5}", "chr", "start", "end", "overlap", "pwd");
        self.blocks.keys().sorted().fold(Ok(()), |_result, chr  | {
            self.blocks[chr].iter().fold(Ok(()), |result, block| {
                result.and_then(|_| writeln!(f, "{}", block))
            })
        })
    }
}


///Chromosome Block for Jackknife resampling. Implemented within struct `JackknifeBlocks`, which is itself implemented
///within structs `pwd_from_stdin::pileup::Comparison`
/// - chromosome  : chromosome name (warn: not by index -> cannot handle specific names, such as "X", "MT", "chr7", etc.)
/// - range       : [start, end[ coordinates of the block.
/// - site_counts : counts the number of overlapping SNPs for a given pair of individuals.
/// - pwd_counts  : counts the number of pairwise differences for a given pair of individuals.
/// 
/// # Traits :
///  - `Eq` and `PartialEq` : with `Self`.
///  - `Eq` and `PartialEq` : with struct `pwd_from_stdin::genome::SNPCoord`
///     - mainly implemented in order to retrieve the block corresponding with a given SNP position.
///       See: `Self::find_block()` or `JackknifeBlocks::find_blocks()`
///  - `Hash`               : along `chromosome` and `range` fields.
///  - `Display`            : Pretty print for file and/or console output Recursively called by `JackknifeBlocks` when it itself is displayed.
#[derive(Debug)]
pub struct JackknifeBlock {
    pub chromosome : u8,
    pub range      : Range<u32>,
    pub site_counts: u32,
    pub pwd_counts : u32,
}

impl JackknifeBlock {
    pub fn new(chromosome: u8, start: u32, end: u32) -> JackknifeBlock {
        JackknifeBlock{chromosome, range: Range{start, end}, site_counts: 0, pwd_counts:0}
    }

    /// Incrementer for the `pwd_counts` field. Called by `pwd_from_stdin::Comparison::compare()`
    /// when a pairwise difference has been declared.
    pub fn add_pwd(&mut self) {
        self.pwd_counts += 1;
    }

    /// Incremented for the `site_counts` field. Infaillibly called by `pwd_from_stdin::Comparison::compare()`
    pub fn add_count(&mut self) {
        self.site_counts += 1;
    }
}

impl fmt::Display for JackknifeBlock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
            "{: <4} - {: <15} - {: <15} - {: <5} - {: <5}",
            self.chromosome,
            self.range.start,
            self.range.end,
            self.site_counts,
            self.pwd_counts
        )
    }
}

impl PartialEq<SNPCoord> for JackknifeBlock {
    fn eq(&self, other: &SNPCoord) -> bool {
        self.chromosome == other.chromosome && self.range.start >= other.position && self.range.end < other.position
    }
    
}

impl PartialEq<JackknifeBlock> for JackknifeBlock {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.range == other.range 
    }
}

impl Eq for JackknifeBlock {}

impl Hash for JackknifeBlock {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.range.hash(state);
    }
}

#[cfg(test)]
mod tests{
    use super::*;
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

}
