use std::{cmp::Ordering, fmt::Display, hash::{Hash, Hasher}};

pub mod chr_index;
pub use chr_index::ChrIdx;

pub mod position;
pub use position::Position;

use crate::jackknife::JackknifeBlock;


pub trait GenomicCoordinate {
    fn coordinate(&self) -> &'_ Coordinate;
}

#[derive(Debug, Clone, Copy)]
pub struct Coordinate {
    pub chromosome: ChrIdx,
    pub position  : Position,
}

impl Display for Coordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&format!("[{: <3} {: >9}]", self.chromosome, self.position), f)
    }
}

impl Coordinate {
    #[must_use]
    pub fn new(chromosome: impl Into<ChrIdx>, position: impl Into<Position>) -> Self {
        Self{chromosome: chromosome.into(), position: position.into()}
    }
}

impl Ord for Coordinate {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.chromosome, self.position).cmp(&(other.chromosome, other.position))
    }
}

impl PartialOrd for Coordinate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq<Self> for Coordinate {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.position == other.position
    }
}

impl PartialEq<JackknifeBlock> for Coordinate {
    fn eq(&self, other: &JackknifeBlock) -> bool {
        other.chromosome == self.chromosome && self.position >= other.range.start && self.position < other.range.end
    }
}

impl Eq for Coordinate {}

impl Hash for Coordinate {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.position.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let (chr, pos) = (12, 230_564_555);
        let want = format!("{:_^18}", format!("[{chr: <3} {pos: >9}]"));
        let got = format!("{:_^18}", Coordinate::new(ChrIdx(chr), Position(pos)));
        println!("want: {want}\ngot : {got}");
        assert_eq!(want, got);
    }
}