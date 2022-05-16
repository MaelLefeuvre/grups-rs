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