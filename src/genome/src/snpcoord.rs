use crate::jackknife::*;
use std::{
    hash::{Hash, Hasher},
    cmp::Ordering,
};

/// A simple struct representing an SNPCoordinate position.
/// Note that :
///   - REF/ALT positions are optional, but are required if the user
///     requires known-sites filtration.
///   - SNPCoordinates implement Equality with JackknifeBlock structs.
///     See: `pwd_from_stdin::jackknife::JackknifeBlock`
///   - Hashable, but only in regards to chr and pos.
///     -> don't insert alternate variations.
#[derive(Debug, Clone)]
pub struct SNPCoord {
    pub chromosome : u8,
    pub position   : u32,
    pub reference  : Option<char>,
    pub alternate  : Option<char>,
}

impl Ord for SNPCoord {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.chromosome, self.position).cmp(&(other.chromosome, other.position))
    }
}

impl PartialOrd for SNPCoord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq<SNPCoord> for SNPCoord {
    fn eq(&self, other: &Self) -> bool { 
        self.chromosome == other.chromosome && self.position == other.position
    }
}

impl PartialEq<JackknifeBlock> for SNPCoord {
    fn eq(&self, other: &JackknifeBlock) -> bool {
        other.chromosome == self.chromosome && self.position >= other.range.start && self.position < other.range.end
    }
}

impl Eq for SNPCoord {}

impl Hash for SNPCoord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chromosome.hash(state);
        self.position.hash(state);
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snp_full_equality() {
        let chromosome1 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        let chromosome2 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }

    #[test]
    fn snp_partial_equality() {
        let chromosome1 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('T'), alternate: Some('G')};
        let chromosome2 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('A'), alternate: Some('C')};
        assert_eq!(chromosome1, chromosome2) 
    }
}
