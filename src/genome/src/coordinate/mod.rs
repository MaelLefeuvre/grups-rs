use std::fmt::Display;

pub mod chr_index;
pub use chr_index::{ChrIdx, ChrIdxError};

pub mod position;
pub use position::{Position, ParsePositionError};

// import for internal use 
extern crate coordinate_derive;
use coordinate_derive::*;

// Derive macro prelude.
pub mod derive {
    extern crate coordinate_derive;
    pub use coordinate_derive::*;
}

/// Padding values for chromosome display
pub const CHR_FORMAT_LEN: usize = 2;
/// Padding values for position display
pub const POS_FORMAT_LEN: usize = 9;

/// Public method trait to access the genomic coordinate of any struct containing a Coordinate struct.
pub trait GenomeCoordinate {
    /// Return a reference to an inner Coordinate struct, representing the genomic location of 
    /// the outer struct
    fn coordinate(&self) -> &'_ Coordinate;

    /// Match the genomic coordinate of the struct with any struct implementing this trait.
    /// # Example
    ///```rust 
    /// # use genome::snp::Allele;
    /// # use genome::coordinate::{Coordinate, ChrIdx, Position, GenomeCoordinate};
    /// # use genome::coordinate::derive::*;
    /// # #[derive(Debug)]
    /// # struct SNP{coordinate: Coordinate, allele: Allele};
    /// # impl GenomeCoordinate for SNP { 
    /// #     fn coordinate(&self) -> &Coordinate { &self.coordinate }
    /// # }
    /// 
    /// let coordinate = Coordinate::new(ChrIdx(10), Position(20));
    /// let my_ref     = SNP{coordinate, allele: Allele::A};
    /// 
    /// assert!(my_ref.matches(&coordinate))
    /// 
    /// ```
    fn matches(&self, other: &impl GenomeCoordinate) -> bool {
        self.coordinate() == other.coordinate()
    }
}

/// Coordinate represents a discrete genomic position, as in 10:591321
/// 
/// Coordinates can be compared and ordered using chromosome and position information. Chromosome, taking priority.
/// # Example: 
/// ```rust
/// use genome::coordinate::{Coordinate, ChrIdx, Position};
/// 
/// let chr_10_200k = Coordinate::new(ChrIdx(10), Position(200_000)); // 10:200000
/// let chr_20_100k = Coordinate::new(ChrIdx(20), Position(100_000)); // 20:100000
/// 
/// assert!(chr_10_200k < chr_20_100k);
/// assert_ne!(chr_10_200k, chr_20_100k);
/// ```
/// Notice the chromosome index takes priority over the position.
/// Thus, the first position is considered _**less**_ than the second, despite having a greater 
/// position (200k vs 100k)
///
/// This behavior can be extended to other structs, using the provided derive macros, as long as the struct houses a Coordinate within `coordinate` feld
/// 
/// # Example: 
/// ```rust 
/// # use genome::snp::Allele;
/// # use genome::coordinate::{Coordinate, ChrIdx, Position, GenomeCoordinate};
/// use genome::coordinate::derive::*;
/// 
/// #[derive(Debug, CoordBorrow, CoordEq)]
/// struct SNP{coordinate: Coordinate, allele: Allele};
///
/// // This trait can also be trivially implemented using derive(Coord)
/// impl GenomeCoordinate for SNP { 
///     fn coordinate(&self) -> &Coordinate { &self.coordinate }
/// }
/// 
/// let coordinate = Coordinate::new(ChrIdx(10), Position(20));
/// let my_ref     = SNP{coordinate, allele: Allele::A};
/// let my_alt     = SNP{coordinate, allele: Allele::T};
/// // Since we derived CoordEq,tThese two value are considered equal, even though their allele differ.
/// assert_eq!(my_ref, my_alt)
/// ```
/// 
#[derive(Debug, Clone, Copy, CoordEq, CoordBlockEq, CoordOrd, CoordHash)]
pub struct Coordinate {
    pub chromosome: ChrIdx,
    pub position  : Position,
}

impl GenomeCoordinate for Coordinate {
    fn coordinate(&self) -> &'_ Coordinate {
        self
    }
    fn matches(&self, other: &impl GenomeCoordinate) -> bool {
        self == other.coordinate()
    }
}

impl Display for Coordinate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&format!("[{: <CHR_FORMAT_LEN$} {: >POS_FORMAT_LEN$}]", self.chromosome, self.position), f)
    }
}

impl Coordinate {
    #[must_use]
    pub fn new(chromosome: impl Into<ChrIdx>, position: impl Into<Position>) -> Self {
        Self{chromosome: chromosome.into(), position: position.into()}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display() {
        let (chr, pos) = (12, 230_564_555);
        let want = format!("{:_^18}", format!("[{chr: <CHR_FORMAT_LEN$} {pos: >POS_FORMAT_LEN$}]"));
        let got  = format!("{:_^18}", Coordinate::new(ChrIdx(chr), Position(pos)));
        println!("want: {want}\ngot : {got}");
        assert_eq!(want, got);
    }
}