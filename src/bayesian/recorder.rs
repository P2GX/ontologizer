use crate::bayesian::state::{CountableState, State};

pub trait Recorder<S: State> {
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, state: &S);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(&mut self);
}

struct Count {
    counts: Vec<u32>,
    n: u32,
}

impl<S> Recorder<S> for Count
where
    S: CountableState,
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

struct Frequency {
    counts: Vec<u32>,
    freqs: Vec<f32>,
    n: u32,
}

impl<S> Recorder<S> for Frequency
where
    S: CountableState,
{
    fn initialize(state: &S) -> Self {
        Self {
            counts: vec![0; state.n_all()],
            freqs: vec![0.0; state.n_all()],
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
        self.freqs = self
            .counts
            .iter()
            .map(|&x| x as f32 / self.n as f32)
            .collect();
    }
}