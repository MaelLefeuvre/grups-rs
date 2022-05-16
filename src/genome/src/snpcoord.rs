use crate::jackknife::JackknifeBlock;
use std::{
    hash::{Hash, Hasher},
    cmp::Ordering,
};

const CHR_FORMAT_LEN: usize = 2;
const POS_FORMAT_LEN: usize = 9;

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

impl std::fmt::Display for SNPCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{: <CHR_FORMAT_LEN$} {: >POS_FORMAT_LEN$}]: {} {}", 
            self.chromosome,
            self.position,
            match self.reference {
                Some(nuc) => nuc.to_string(),
                None => "None".to_string(),
            },
            match self.alternate {
                Some(nuc) => nuc.to_string(),
                None => "None".to_string(),
            },        
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{self, prelude::SliceRandom};

    const N_ITERS: u32 = 1_000_000;

    #[test]
    fn snpcoord_ordering(){
        let (chromosome, position) = (10, 100_000);
        let coord = SNPCoord {chromosome, position, reference: Some('A'), alternate: Some('C')};
        let deviations: [i32; 3] = [-1, 0, 1];

        for nucleotides in [('A', 'T'), ('C', 'G')] { // Differing reference and/or alternate nucleotides should not impact ordering.
            let (reference, alternate) = (Some(nucleotides.0), Some(nucleotides.1));
            for chromosome_deviation in deviations {
                for position_deviation in deviations {
                    let chromosome = chromosome as i32 + chromosome_deviation;
                    let position   = position   as i32 + position_deviation;

                    let other_coord = SNPCoord {chromosome: chromosome as u8, position: position as u32, reference, alternate};

                    match chromosome_deviation {
                        -1 => assert!(other_coord < coord),
                         1 => assert!(other_coord > coord),
                         0 => match position_deviation {
                            -1 => assert!(other_coord < coord),
                             1 => assert!(other_coord > coord),
                             0 => assert_eq!(other_coord, coord),
                             _ => continue
                         }
                         _ => continue
                    }
                }
            }
        }

    }

    #[test]
    fn snpcoord_full_equality() {
        let coord1 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        let coord2 = SNPCoord {chromosome: 1, position: 100510, reference: Some('A'), alternate: Some('C')};
        assert_eq!(coord1, coord2) 
    }

    #[test]
    fn snpcoord_partial_equality() {
        let coord1 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('T'), alternate: Some('G')};
        let coord2 = SNPCoord {chromosome: 2, position: 16541561, reference: Some('A'), alternate: Some('C')};
        assert_eq!(coord1, coord2) 
    }

    #[test]
    fn jackknife_equality() {
        let coord= SNPCoord{chromosome: 10, position: 11_000, reference: None, alternate: None};

        //snpcoord == block if block.start <= snpcoord.position < block.end AND block.chr == snpcoord.chr
        assert_eq!(coord, JackknifeBlock::new(10, 11_000, 12_000)); // le. range    ; same chromosome.
        assert_eq!(coord, JackknifeBlock::new(10, 10_000, 12_000)); // within range ; same chromosome.
    }

    #[test]
    fn jackknife_inequality() {
        let coord= SNPCoord{chromosome: 10, position: 11_000, reference: None, alternate: None};

        assert_ne!(coord, JackknifeBlock::new(10, 11_500, 12_000)); // lt. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(10, 10_000, 11_000)); // ge. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(10, 10_000, 10_500)); // gt. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(11, 10_000, 12_000)); // within range ; different chromosome.
    }

    #[test]
    fn hash_block() {
        let mut rng = rand::thread_rng();

        let nucleotides = vec!['A', 'C', 'G', 'T'];
        let mut test_hashset = std::collections::HashSet::new();
        for chromosome in 1..22 {
            for position in (1..N_ITERS).step_by(1000) {
                let mut random_nucl = nucleotides.choose_multiple(&mut rng, 2).map(|nuc| *nuc);
                let coord = SNPCoord{chromosome, position, reference: random_nucl.next(), alternate: random_nucl.next()};
                assert!(test_hashset.insert(coord));
            }
        }

        for chromosome in 1..22 {
            for position in (1..N_ITERS).step_by(1000) {
                let coord = SNPCoord{chromosome, position, reference: None, alternate: None};
                let retrieved_coord = test_hashset.get(&coord);
                assert!(retrieved_coord.is_some());
            }
        }
    }


    #[test]
    fn display_some_refalt() {
        let (chromosome, position) = (1, 1_357_165);

        let coord = SNPCoord {chromosome, position, reference: Some('A'), alternate: Some('C')};
        let expected_output = format!("[{: <CHR_FORMAT_LEN$} {: >POS_FORMAT_LEN$}]: {} {}", chromosome, position, 'A', 'C');
        assert_eq!(format!("{coord}"), expected_output)
    }

    #[test]
    fn display_none_refalt() {
        let (chromosome, position) = (1, 1_357_165);

        let coord = SNPCoord {chromosome, position, reference: None, alternate: None};
        let expected_output = format!("[{: <CHR_FORMAT_LEN$} {: >POS_FORMAT_LEN$}]: {} {}", chromosome, position, "None", "None");
        assert_eq!(format!("{coord}"), expected_output)
    }

}
