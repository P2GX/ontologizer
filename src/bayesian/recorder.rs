use crate::bayesian::state::{BinaryParameterState, State};
use crate::core::result::Measure;

pub trait Recorder<S: State> {
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, state: &S);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(&mut self);
}

pub struct Count {
    counts: Vec<u32>,
    n: u32,
}

impl Count {
    pub fn iter(&self) -> std::slice::Iter<u32> {
        self.counts.iter()
    }

    pub fn n(&self) -> u32 {
        self.n
    }
}

impl<S> Recorder<S> for Count
where
    S: BinaryParameterState,
{
    fn initialize(state: &S) -> Self {
        Self {
            counts: vec![0; state.n_all()],
            n: 0,
        }
    }
    fn record(&mut self, state: &S) {
        for i in 0..state.n_all() {
            if state.get(i) {
                self.counts[i] += 1;
            }
        }
        self.n += 1;
    }

    fn finalize(&mut self) {}
}

impl Measure for Count {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.counts.iter().map(|&c| c as f64)
    }
}

pub struct Frequency {
    counts: Vec<u32>,
    frequencies: Vec<f32>,
    n: u32,
}

impl<S> Recorder<S> for Frequency
where
    S: BinaryParameterState,
{
    fn initialize(state: &S) -> Self {
        Self {
            counts: vec![0; state.n_all()],
            frequencies: vec![0.0; state.n_all()],
            n: 0,
        }
    }
    fn record(&mut self, state: &S) {
        for i in 0..state.n_all() {
            if state.get(i) {
                self.counts[i] += 1;
            }
        }
        self.n += 1;
    }

    fn finalize(&mut self) {
        self.frequencies = self
            .counts
            .iter()
            .map(|&x| x as f32 / self.n as f32)
            .collect();
    }
}

// ... existing struct Frequency ...

impl Measure for Frequency {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.frequencies.iter().map(|&f| f as f64)
    }
}
