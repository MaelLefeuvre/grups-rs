use rand::Rng;
use std::fmt::Debug;
use std::ops::Range;
use rand::distributions::uniform::SampleUniform;
use std::cmp::PartialOrd;

/// Trait defining a pedigree parameter. This struct is mainly leveraged by `super::ParamRateGenerator` to generate
/// constant/random values, according to the user-provided rates.
pub trait PedParam<T>: std::fmt::Debug {
    fn value(&mut self) -> T;
}

impl<T: 'static + Copy + SampleUniform + PartialOrd + Debug> dyn PedParam<T> {
    /// Instantiate a new PedParam, given a vector of values.
    /// 
    /// # Arguments:
    /// - `vec`: vector of user-provided value(s). Must either be of len 1 or 2.
    /// 
    /// # Panics
    /// - whenever `vec.len() > 2`

    pub fn from_vec(vec: &[T]) -> Box<dyn PedParam<T>> {
        match vec.len() {
            // If length is 1, the user is requesting a set, constant rate
            // => instantiate a `PedParamConst`
            1 => Box::new(PedParamConst::new(vec[0])),                    

            // If length is 2, the user is requesting for random values within a range 
            // => instantiate a `PedParamRange`
            2 => Box::new(PedParamRange::new(vec[0], vec[1])),

            // Uh oh... This should never happen, as sanity checks are performed upstream. 
            _ => panic!("Failed to parse pedigree parameter!"),
        }
    }
}

/// Constant pedigree parameter. Calling `self.value()` on this struct will always return the same value (the one given by the user)
/// # Fields:
/// - `inner` : user-provided rate.
#[derive(Debug)]
struct PedParamConst<T> {
    inner: T
}

impl<T: Copy + Debug> PedParam<T> for PedParamConst<T> {
    /// Return the rate previously given by the user.
    fn value(&mut self) -> T {
        self.inner
    }
}

impl<T> PedParamConst<T> {
    /// Instantiate a new constant PedParam. 
    pub fn new(value: T) -> Self {
        PedParamConst{inner: value}
    }
}


/// Random pedigree parameter. Calling `self.value()` on this struct will return a random value within a user-provided range.
/// # Fields:
/// - `inner`: random number generator.
/// - `range`: user-defined range of acceptable rates.
#[derive(Debug)]
struct PedParamRange<T> {
    inner: rand::prelude::ThreadRng,
    range: Range<T>
}

impl<T: SampleUniform + PartialOrd + Copy + Debug> PedParam<T> for PedParamRange<T> {
    /// Return a randomly generated rate, within the user-provided range.
    fn value(&mut self) -> T {
        self.inner.gen_range(self.range.start..=self.range.end)
    }
}

impl<T> PedParamRange<T> {
    /// Instantiate a new random PedParam.
    pub fn new(start: T, end: T) -> Self {
        let inner = rand::thread_rng();
        let range = Range{start, end};
        Self {inner, range}
    }
}