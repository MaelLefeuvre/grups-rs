use rand::distributions::uniform::SampleUniform;
use std::cmp::PartialOrd;
use std::fmt::{Debug, Display};

use super::PedParam;

/// Wrapper struct containing two PedParam traits. `self.inner[i]` = PedParam of pair[i]
/// 
/// Each `PedParam` can bear two different behaviors, according to the type of user-provided data
/// - if the provided data is a single value (e.g.: '5' ), `PedParam` is of type `PedParamConst`
///   -> `ParamRateGenerator` will always return a constant value for this `PedParam`(the one provided by the user)
///      when calling `gen_random_values()`
/// - if the provided data is a range (e.g.: '5-10'), `PedParam` is of type `PedParamRange`
///   -> `ParamRateGenerator` will return a random value for this `PedParam` (within the provided value)
///      when calling | gen_random_values()`
#[derive(Debug)]
pub struct ParamRateGenerator<T>{
    inner: [Box<dyn PedParam<T>>; 2],
}

impl<T: Display> Display for ParamRateGenerator<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} vs. {}", self.inner[0], self.inner[1])
    }
}

impl<T: Display> ParamRateGenerator<T> {
    /// Instantiate a new generartor from a given user-input.
    /// # Arguments:
    /// - `rates`   : vector of user provided rate (e.g. contamination rate, sequencing error rate, etc.)
    ///               each inner vector must either be of len 1 or 2. any value > 2 will cause a panic.
    /// - `indices` : pileup indices of the two individuals being compared. (See `pwd_from_stdin::Comparison::get_pair_indices()`)
    /// 
    /// # Panics:
    /// - whenever `rates[i].len()` > 2
    /// - if self.inner.len() > 2
    pub fn from_user_input(rates: &[Vec<T>], indices: [usize; 2]) -> Self 
    where 
        T: 'static + Copy + SampleUniform + PartialOrd + Debug
    {
        let mut inner = Vec::new();
        // ---- Instantiate two PedParam. One for each compaired individual. 
        for index in &indices {
            let rate_idx = index % rates.len() ; // Wrap around values if the user did not provide enough ranges!
            let rate_gen: Box<dyn PedParam<T>> = <dyn PedParam<T>>::from_vec(&rates[rate_idx]);
            inner.push(rate_gen);
        }
        let inner: [Box<dyn PedParam<T>>; 2 ] = inner.try_into()
            .expect("Invalid ParamRateGenerator array length.");
            
        Self {inner}
    }

    /// Generate two random values, one for each compaired individual.
    pub fn gen_random_values(&mut self) -> [T; 2] {
        [self.inner[0].value(), self.inner[1].value()]
    }
}
