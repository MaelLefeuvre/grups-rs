use rand::distributions::uniform::SampleUniform;
use std::cmp::PartialOrd;

use super::PedParam;

pub struct ParamRateGenerator<T>{
    inner: [Box<dyn PedParam<T>>; 2],
}

impl<T> ParamRateGenerator<T> {
    pub fn from_user_input(rates: &Vec<Vec<T>>, indices: [usize; 2]) -> Self 
    where 
        T: 'static + Copy + SampleUniform + PartialOrd
    {
        let mut inner = Vec::new();
        for i in 0..2 {
            let rate_idx = indices[i] % rates.len() ; // Wrap around values if the user did not provide enough ranges!
            let rate_gen: Box<dyn PedParam<T>> = <dyn PedParam<T>>::from_vec(&rates[rate_idx]);
            inner.push(rate_gen);
        }
        let inner: [Box<dyn PedParam<T>>; 2 ] = inner.try_into().unwrap_or_else(|v: Vec<Box<dyn PedParam<T>>> | panic!("Expected a Vec of length {} but it was {}", 2, v.len()));
        Self {inner}
    }

    pub fn gen_random_values(&mut self) -> [T; 2] {
        [self.inner[0].value(), self.inner[1].value()]
    }
}
