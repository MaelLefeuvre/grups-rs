use rand::Rng;
use rand::RngCore;
use std::ops::Range;
use rand::distributions::uniform::SampleUniform;
use std::cmp::PartialOrd;

pub trait PedParam<T>{
    fn value(&mut self) -> T;
}

impl<T: 'static + Copy + SampleUniform + PartialOrd> dyn PedParam<T> {
    pub fn from_vec(vec: &Vec<T>) -> Box<dyn PedParam<T>> {
        match vec.len() {
            1 => Box::new(PedParamConst::new(vec[0])),
            2 => Box::new(PedParamRange::new(vec[0], vec[1])),
            _ => panic!("Failed to parse pedigree parameter!"),
        }
    }
}

struct PedParamConst<T> {
    inner: T
}

impl<T: Copy> PedParam<T> for PedParamConst<T> {
    fn value(&mut self) -> T {
        self.inner
    }
}

impl<T> PedParamConst<T> {
    pub fn new(value: T) -> Self {
        PedParamConst{inner: value}
    }
}

struct PedParamRange<T> {
    inner: Box<dyn RngCore>,
    range: Range<T>
}

impl<T: SampleUniform + PartialOrd + Copy> PedParam<T> for PedParamRange<T> {
    fn value(&mut self) -> T {
        self.inner.gen_range(self.range.start, self.range.end)
    }
}

impl<T> PedParamRange<T> {
    pub fn new(start: T, end: T) -> Self {
        let inner = Box::new(rand::thread_rng());
        let range = Range{start, end};
        Self {inner, range}
    }
}

#[derive(Debug, Clone)]
pub struct PedigreeParams {
    pub snp_downsampling_rate : f64,
    pub af_downsampling_rate  : f64,
    pub seq_error_rate        : [f64; 2],
    pub contam_rate           : [f64; 2],
}

impl PedigreeParams {
    pub fn new(snp_downsampling_rate: f64, af_downsampling_rate: f64, seq_error_rate: [f64; 2], contam_rate: [f64; 2]) -> Self {
        PedigreeParams{snp_downsampling_rate, af_downsampling_rate, seq_error_rate, contam_rate}
    }
}

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

