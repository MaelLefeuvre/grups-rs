mod allele;
pub use allele::Allele;
pub use allele::ParseAlleleError;


use anyhow::Result;
use located_error::LocatedError;
use std::error::Error;

use crate::coordinate::{Coordinate, ChrIdx, Position, derive::*};

//pub use crate::snp::Allele;

/// A simple struct representing an `SNPCoordinate` position.
/// Note that :
///   - REF/ALT positions are optional, but are required if the user
///     requires known-sites filtration.
///   - `SNPCoordinates` implement Equality with `JackknifeBlock` structs.
///     See: `pwd_from_stdin::jackknife::JackknifeBlock`
///   - Hashable, but only in regards to chr and pos.
///     -> don't insert alternate variations.
#[derive(Debug, Clone, Copy, Coord, CoordEq, CoordBlockEq, CoordOrd, CoordHash, CoordBorrow)]
pub struct SNPCoord {
    pub coordinate: Coordinate,
    pub reference  : Allele,
    pub alternate  : Allele,
}

impl std::fmt::Display for SNPCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {} {}", 
            self.coordinate(),
            self.reference,
            self.alternate,
        )
    }
}

impl SNPCoord {

    pub fn new(chromosome: impl Into<ChrIdx>, position: impl Into<Position>, reference: Allele, alternate: Allele) -> Self {
        let coordinate = Coordinate::new(chromosome, position);
        Self{coordinate, reference, alternate}
    }

    pub fn try_new<T>(chromosome: impl Into<ChrIdx>, position: impl Into<Position>, reference: T, alternate: T) -> Result<Self> 
    where   T: TryInto<Allele>,
            T::Error: Error + Sync + Send + 'static
    {
        let coordinate = Coordinate::new(chromosome, position);
        let reference = reference.try_into().with_loc(|| format!("While parsing coordinate {coordinate}"))?;
        let alternate = alternate.try_into().with_loc(|| format!("While parsing coordinate {coordinate}"))?;
        Ok(Self{coordinate, reference, alternate})
    }
    
    /// Check if this coordinate contains known REF/ALT alleles.
    #[must_use]
    pub fn has_known_alleles(&self) -> bool {
        self.reference.is_known() && self.alternate.is_known()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use  crate::jackknife::JackknifeBlock;
    use rand::{self, prelude::SliceRandom};
    use anyhow::Result;
    const N_ITERS: u32 = 1_000_000;

    #[test]
    fn snpcoord_ordering() -> Result<()> {
        let (chromosome, position) = (10, 100_000);
        let coord = SNPCoord::try_new(chromosome, position, 'A', 'C')?;
        let deviations: [i32; 3] = [-1, 0, 1];

        for nucleotides in [('A', 'T'), ('C', 'G')] { // Differing reference and/or alternate nucleotides should not impact ordering.
            let (reference, alternate) = (nucleotides.0, nucleotides.1);
            for chromosome_deviation in deviations {
                for position_deviation in deviations {
                    let chromosome = chromosome as i32 + chromosome_deviation;
                    let position   = position   as i32 + position_deviation;

                    let other_coord = SNPCoord::try_new(chromosome as u8, position as u32, reference, alternate)?;

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
        Ok(())
    }

    #[test]
    fn snpcoord_full_equality() -> Result<()> {
        let coord1 = SNPCoord::try_new(1, 100510, 'A', 'C')?;
        let coord2 = SNPCoord::try_new(1, 100510, 'A', 'C')?;
        assert_eq!(coord1, coord2);
        Ok(())
    }

    #[test]
    fn snpcoord_partial_equality() -> Result<()> {
        let coord1 = SNPCoord::try_new(2, 16541561, 'T', 'G')?;
        let coord2 = SNPCoord::try_new(2, 16541561, 'A', 'C')?;
        assert_eq!(coord1, coord2);
        Ok(())
    }

    #[test]
    fn jackknife_equality() -> Result<()> {
        let coord = SNPCoord::try_new(10, 11_000, 'N', 'N')?;

        //snpcoord == block if block.start <= snpcoord.position < block.end AND block.chr == snpcoord.chr
        assert_eq!(coord, JackknifeBlock::new(10, 11_000, 12_000)); // le. range    ; same chromosome.
        assert_eq!(coord, JackknifeBlock::new(10, 10_000, 12_000)); // within range ; same chromosome.
        Ok(())
    }

    #[test]
    fn jackknife_inequality() -> Result<()>  {
        let coord = SNPCoord::try_new(10, 11_000, 'N', 'N')?;

        assert_ne!(coord, JackknifeBlock::new(10, 11_500, 12_000)); // lt. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(10, 10_000, 11_000)); // ge. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(10, 10_000, 10_500)); // gt. range    ; same chromosome
        assert_ne!(coord, JackknifeBlock::new(11, 10_000, 12_000)); // within range ; different chromosome.
        Ok(())
    }

    #[test]
    fn hash_block() -> Result<()> {
        let mut rng = rand::thread_rng();

        let nucleotides = vec!['A', 'C', 'G', 'T'];
        let mut test_hashset = std::collections::HashSet::new();
        for chromosome in 1..22 {
            for position in (1..N_ITERS).step_by(1000) {
                let mut random_nucl = nucleotides.choose_multiple(&mut rng, 2);
                let coord = SNPCoord::try_new(chromosome, position, *random_nucl.next().unwrap(), *random_nucl.next().unwrap())?;
                assert!(test_hashset.insert(coord));
            }
        }

        for chromosome in 1..22 {
            for position in (1..N_ITERS).step_by(1000) {
                let coord = SNPCoord::try_new(chromosome, position, 'N', 'N')?;
                let retrieved_coord = test_hashset.get(&coord);
                assert!(retrieved_coord.is_some());
            }
        }
        Ok(())
    }


    #[test]
    fn display_some_refalt() -> Result<()> {
        use crate::coordinate::{CHR_FORMAT_LEN, POS_FORMAT_LEN};

        let (chromosome, position) = (1, 1_357_165);
        let coord = SNPCoord::try_new(chromosome, position, 'A', 'C')?;
        let expected_output = format!("[{chromosome: <CHR_FORMAT_LEN$} {position: >POS_FORMAT_LEN$}]: A C");
        assert_eq!(format!("{coord}"), expected_output);
        Ok(())
    }

    #[test]
    fn display_none_refalt() -> Result<()> {
        use crate::coordinate::{CHR_FORMAT_LEN, POS_FORMAT_LEN};
        let (chromosome, position) = (1, 1_357_165);

        let coord = SNPCoord::try_new(chromosome, position, 'N', 'N')?;
        let expected_output = format!("[{chromosome: <CHR_FORMAT_LEN$} {position: >POS_FORMAT_LEN$}]: N N");
        assert_eq!(format!("{coord}"), expected_output);
        Ok(())
    }

}
