use coordinate_derive::{Coord, CoordBorrow, CoordEq, CoordOrd, CoordHash};
use genome::coordinate::{Coordinate};
use genome::coordinate::{ChrIdx, Position};

#[derive(Debug, Coord, CoordEq, CoordHash, CoordBorrow, CoordOrd)]
struct TestMacro{coordinate: Coordinate}

impl TestMacro {
    fn new(chromosome: ChrIdx, position: Position) -> Self {
        Self{coordinate: Coordinate{chromosome, position}}
    }
}

#[test]
fn create_derive() {
    #[derive(Coord, CoordEq, CoordOrd, CoordHash, CoordBorrow)]
    struct TestMacro{coordinate: Coordinate}

    
}

#[test]
fn derive_equal() {
    let x = TestMacro::new(ChrIdx(10), Position(250_000));
    let y = TestMacro::new(ChrIdx(11), Position(150_000));
    assert_eq!(x, x);
    assert_ne!(x, y);

    assert!(y > x);
}

#[test]
fn derive_ord() {
    let chr_10_150 = TestMacro::new(ChrIdx(10), Position(150_000));
    let chr_10_250 = TestMacro::new(ChrIdx(10), Position(250_000));
    let chr_11_150 = TestMacro::new(ChrIdx(11), Position(150_000));
    let chr_11_250 = TestMacro::new(ChrIdx(11), Position(250_000));

    // ---- position 10:150
    // vs. 10:150: 10 = 10 and 150 < 250
    assert!(!chr_10_150.lt(&chr_10_150)); // False
    assert!( chr_10_150.le(&chr_10_150)); // True 
    assert!( chr_10_150.ge(&chr_10_150)); // True 
    assert!(!chr_10_150.gt(&chr_10_150)); // False
    // vs. 10:250: 10 = 10 and 150 < 250
    assert!( chr_10_150.lt(&chr_10_250)); // True  
    assert!( chr_10_150.le(&chr_10_250)); // True   
    assert!(!chr_10_150.ge(&chr_10_250)); // False 
    assert!(!chr_10_150.gt(&chr_10_250)); // False 
    // vs. 11:150: 10 < 11                                      
    assert!( chr_10_150.lt(&chr_11_150)); // True 
    assert!( chr_10_150.le(&chr_11_150)); // True 
    assert!(!chr_10_150.ge(&chr_11_150)); // False
    assert!(!chr_10_150.gt(&chr_11_150)); // False
    // vs. 11:250: 10 < 11                                      
    assert!( chr_10_150.lt(&chr_11_250)); // True 
    assert!( chr_10_150.le(&chr_11_250)); // True 
    assert!(!chr_10_150.ge(&chr_11_250)); // False
    assert!(!chr_10_150.gt(&chr_11_250)); // False

    // ---- position 10:250
    // vs. 10:250: 10 = 10 and 250 = 250
    assert!(!chr_10_250.lt(&chr_10_250)); // False
    assert!( chr_10_250.le(&chr_10_250)); // True 
    assert!( chr_10_250.ge(&chr_10_250)); // True 
    assert!(!chr_10_250.gt(&chr_10_250)); // False
    // vs. 11:150: 10 < 11
    assert!( chr_10_250.lt(&chr_11_150)); // True 
    assert!( chr_10_250.le(&chr_11_150)); // True 
    assert!(!chr_10_250.ge(&chr_11_150)); // False
    assert!(!chr_10_250.gt(&chr_11_150)); // False
    // vs. 11:250: 10 < 11
    assert!( chr_10_250.lt(&chr_11_250)); // True 
    assert!( chr_10_250.le(&chr_11_250)); // True 
    assert!(!chr_10_250.ge(&chr_11_250)); // False
    assert!(!chr_10_250.gt(&chr_11_250)); // False

    // ---- position 11:150
    // vs 11:150: 11 = 11 and 150 = 150
    assert!(!chr_11_150.lt(&chr_11_150)); // False 
    assert!( chr_11_150.le(&chr_11_150)); // True 
    assert!( chr_11_150.ge(&chr_11_150)); // True
    assert!(!chr_11_150.gt(&chr_11_150)); // False
    // vs 11::250: 150 < 250
    assert!( chr_11_150.lt(&chr_11_250)); // True 
    assert!( chr_11_150.le(&chr_11_250)); // True 
    assert!(!chr_11_150.ge(&chr_11_250)); // False
    assert!(!chr_11_150.gt(&chr_11_250)); // False

    // ---- position 11:250
    // vs 11:150: 11 = 11 and 150 = 150
    assert!(!chr_11_250.lt(&chr_11_250)); // True 
    assert!( chr_11_250.le(&chr_11_250)); // True 
    assert!( chr_11_250.ge(&chr_11_250)); // True
    assert!(!chr_11_250.gt(&chr_11_250)); // False
}
