use std::{
    hash::{Hash, Hasher},
    ops::{Range},
    borrow::Borrow
};

// ---- Constant values
const ONE_MORGAN  : f64 = 100.0; 
const ONE_MEGABASE: f64 = 1_000_000.0;

/// Genetic Recombination probabily over a given range.
/// # Fields
/// - `range`: 0-based <start> and <end> coordinates on wich `self.prob` is tied.
/// - `prob` : probability of genetic recombination 
#[derive(Debug, Clone)]
pub struct RecombinationRange {
    range: Range<u32>,
    prob : f64,
}

impl RecombinationRange {
    /// Instantiate a new `RecombinationRange` and probability from a genetic recombination map entry.
    /// # Arguments
    /// - `start`: 0-based start-position of the range. (field 0 of a genetic map file)
    /// - `end`  : 0-based end   position of the range. (field 1 of a genetic map file)
    /// - `rate` : recombination rate observed within [start, end[ (field 2 of a genetic map file)
    pub fn new(start: u32, end: u32, rate: f64) -> RecombinationRange {
        let prob: f64 = rate / ONE_MORGAN / ONE_MEGABASE;  // rate/cM/Mb.
        RecombinationRange{range: Range{start, end}, prob}
    }

    /// Return the probability of a recombination occuring within this range.
    pub fn prob(&self) -> &f64 {
        &self.prob
    }

    /// Convert `self.prob` back to a recombination rate.
    #[cfg(test)]
    pub fn rate(&self) -> f64 {
        self.prob * ONE_MORGAN * ONE_MEGABASE
    }
}

impl PartialEq<RecombinationRange> for RecombinationRange {
    fn eq(&self, other: &Self) -> bool { 
        self.range == other.range
    }
}

impl PartialEq<Range<u32>> for RecombinationRange {
    fn eq(&self, other: &Range<u32>) -> bool {
        self.range == *other
    }
}

impl Eq for RecombinationRange {}

impl Hash for RecombinationRange {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.range.hash(state);
    }
}

impl Borrow<Range<u32>> for RecombinationRange {
    fn borrow(&self) -> &Range<u32> {
        &self.range
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use float_cmp::approx_eq;
    const ULPS: i64 = 2;

    const N_ITERS: u32 = 1_000_000;

    #[test]
    fn prob() {
        for step in 0..N_ITERS {
            let test_rate = f64::from(step) / f64::from(N_ITERS) ;
            let test_prob = test_rate / ONE_MORGAN / ONE_MEGABASE;
            let range = RecombinationRange::new(0, 1000, test_rate);
            assert!( approx_eq!(f64, *range.prob(), test_prob, ulps = ULPS) );
        }
    }

    #[test]
    fn prob_to_rate() {
        for step in 0..N_ITERS {
            let test_rate = f64::from(step) / f64::from(N_ITERS) ;
            let range = RecombinationRange::new(0, 1000, test_rate);
            assert!( approx_eq!(f64, range.rate(), test_rate, ulps = ULPS) );
        }
    }

    #[test]
    fn recomb_equality() {
        for step in 0..N_ITERS {
            let test_rate = f64::from(step) / f64::from(N_ITERS) ;
            let range_1 = RecombinationRange::new(0, 1000, test_rate);
            let range_2 = RecombinationRange::new(0, 1000, test_rate);
            assert_eq!(range_1, range_2);
        }
    }

    #[test]
    fn recomb_partial_equality() {
        // Recomb range should still be equal even if test rates are different.
        for step in 0..N_ITERS  {
            let test_rate_1 = f64::from(step) / f64::from(N_ITERS);
            let range_1 = RecombinationRange::new(0, 1000, test_rate_1);

            let test_rate_2 = test_rate_1 + (1.0 / f64::from(N_ITERS));
            let range_2 = RecombinationRange::new(0, 1000, test_rate_2);
            assert_eq!(range_1, range_2);
        }
    }

    #[test]
    fn recomb_inequality() {
        let range_1 = RecombinationRange::new(0, 1000, 0.5);
        let range_2 = RecombinationRange::new(0, 1001, 0.5);
        assert!(range_1 != range_2);
    }

    #[test]
    fn range_equality() {
        let ranges: Vec<u32> = (0..N_ITERS).step_by(1000).collect();
        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);
            let recomb_range = RecombinationRange::new(start, end, 0.5);
            let range = Range{start, end};
            assert_eq!(recomb_range, range);
        }
    }

    #[test]
    fn range_inequality() {
        let ranges: Vec<u32> = (1..=N_ITERS).step_by(1000).collect();
        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);
            let recomb_range = RecombinationRange::new(start, end, 0.5);
            assert!(recomb_range != Range{start, end: end+1});
            assert!(recomb_range != Range{start, end: end-1});

            assert!(recomb_range != Range{start:start+1, end});
            assert!(recomb_range != Range{start:start-1, end});

            assert!(recomb_range != Range{start:start+1, end: end+1});
            assert!(recomb_range != Range{start:start-1, end: end-1});
        }
    }

    #[test]
    fn borrow_range() {
        let mut test_hashset = HashSet::new();
        let ranges: Vec<u32> = (0..N_ITERS).step_by(1000).collect();
        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);
            let range = RecombinationRange::new(start, end, 0.5);
            assert!(test_hashset.insert(range));
        }

        for step in ranges.windows(2)  {
            let (start, end) = (step[0], step[1]);
            assert!(test_hashset.contains(&Range{start, end}));
        }
    }


}