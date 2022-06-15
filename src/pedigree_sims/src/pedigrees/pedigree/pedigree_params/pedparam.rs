use rand::Rng;
use rand::RngCore;
use std::ops::Range;
use rand::distributions::uniform::SampleUniform;
use std::cmp::PartialOrd;

/// Trait defining a pedigree parameter. This struct is mainly leveraged by `super::ParamRateGenerator` to generate
/// constant/random values, according to the user-provided rates.
pub trait PedParam<T>{
    fn value(&mut self) -> T;
}

impl<T: 'static + Copy + SampleUniform + PartialOrd> dyn PedParam<T> {
    /// Instantiate a new PedParam, given a vector of values.
    /// 
    /// # Arguments:
    /// - `vec`: vector of user-provided value(s). Must either be of len 1 or 2.
    /// 
    /// # Panics
    /// - whenever `vec.len() > 2`

    pub fn from_vec(vec: &Vec<T>) -> Box<dyn PedParam<T>> {
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
struct PedParamConst<T> {
    inner: T
}

impl<T: Copy> PedParam<T> for PedParamConst<T> {
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
struct PedParamRange<T> {
    inner: Box<dyn RngCore>,
    range: Range<T>
}

impl<T: SampleUniform + PartialOrd + Copy> PedParam<T> for PedParamRange<T> {
    /// Return a randomly generated rate, within the user-provided range.
    fn value(&mut self) -> T {
        self.inner.gen_range(self.range.start, self.range.end)
    }
}

impl<T> PedParamRange<T> {
    /// Instantiate a new random PedParam.
    pub fn new(start: T, end: T) -> Self {
        let inner = Box::new(rand::thread_rng());
        let range = Range{start, end};
        Self {inner, range}
    }
}