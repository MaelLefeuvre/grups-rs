use crate::genome::*;
use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::ops::Range;
use itertools::Itertools;

#[derive(Debug)]
pub struct JackknifeBlocks {
    blocks : HashMap<u8, Vec<JackknifeBlock>>
}

impl JackknifeBlocks {
    pub fn new(genome: &Vec<Chromosome>, blocksize: u32) -> JackknifeBlocks {
        let mut jackknives = HashMap::new();
        for chr in genome {
            let mut blocks = Vec::new();
            let mut add_block = |chr: u8, start: u32, end: u32| {
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

    pub fn find_block(&mut self, coordinate: &SNPCoord) -> &mut JackknifeBlock {
        self.blocks
            .get_mut(&coordinate.chromosome)
            .unwrap()
            .iter_mut()
            .find(|block| block.range.contains(&coordinate.position)).unwrap()
    }

    pub fn print(&self) -> () {
        println!("{: <4} - {: <15} - {: <15} - {: <5} - {: <5}", "chr", "start", "end", "overlap", "pwd");
        for chr in self.blocks.keys().sorted() {
            for block in self.blocks[chr].iter(){
                println!("{: <4} - {: <15} - {: <15} - {: <5} - {: <5}",
                         block.chromosome,
                         block.range.start,
                         block.range.end,
                         block.site_counts,
                         block.pwd_counts
                );
            }
        }
        
    }
}

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

    pub fn add_pwd(&mut self) {
        self.pwd_counts += 1;
    }

    pub fn add_count(&mut self) {
        self.site_counts += 1;
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

    #[test]
    fn jackknife_blocks_check_len_unequal() {
       let genome = vec![Chromosome{index:0, name: 1, length: 249250621},];
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[0].name].len(),249251);
    }
    
    #[test]
    fn jackknife_blocks_check_len_equal () {
       let genome = vec![Chromosome{index:0, name: 1, length: 2000001},];
       let blocks = JackknifeBlocks::new(&genome, 1000);
       assert_eq!(blocks.blocks[&genome[0].name].len(), 2000);
    }

}
