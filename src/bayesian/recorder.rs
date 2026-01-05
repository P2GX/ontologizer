use crate::bayesian::proposer::ToggleSwap;
use crate::bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::bayesian::state::{BinaryParameterState, State};
use crate::core::result::Measure;

pub trait Recorder<S: State> {
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, m: &<S as State>::Move, step: usize);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(&mut self, final_step: usize);
}

pub struct Count {
    counts: Vec<usize>,
    swaps: Vec<usize>,
    active_since: Vec<Option<usize>>,
}

impl Count {
    pub fn iter(&self) -> std::slice::Iter<usize> {
        self.counts.iter()
    }

    fn toggle_term(&mut self, idx: usize, step: usize) {
        self.swaps[idx] += 1;
        match self.active_since[idx] {
            Some(start) => {
                // Was ON, turning OFF. Record the duration.
                self.counts[idx] += step - start;
                self.active_since[idx] = None;
            }
            None => {
                // Was OFF, turning ON. Start the timer.
                self.active_since[idx] = Some(step);
            }
        }
    }
}

impl<S> Recorder<S> for Count
where
    S: BinaryParameterState<Move = ToggleSwap>,
{
    fn initialize(state: &S) -> Self {
        let n = state.n_all();
        let mut active_since = vec![None; n];

        for i in 0..n {
            if state.get(i) {
                active_since[i] = Some(0);
            }
        }

        Self {
            counts: vec![0; n],
            swaps: vec![0; n],
            active_since,
        }
    }

    fn record(&mut self, m: &S::Move, step: usize) {
        match *m {
            Toggle(i) => self.toggle_term(i, step),
            Swap(i, j) => {
                self.toggle_term(i, step);
                self.toggle_term(j, step);
            }
        }
    }

    fn finalize(&mut self, final_step: usize) {
        for (i, start_opt) in self.active_since.iter().enumerate() {
            if let Some(start) = start_opt {
                self.counts[i] += final_step - start;
            }
        }
    }
}

impl Measure for Count {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.counts.iter().map(|&x| x as f64)
    }

    fn diagnostics(&self) -> impl Iterator<Item = Option<String>> {
        self.swaps.iter().map(|&x| Some(x.to_string()))
    }

    fn get_score(&self, i: usize) -> f64 {
        self.counts[i] as f64
    }

    fn get_diagnostics(&self, i: usize) -> Option<f64> {
        Some(self.swaps[i] as f64)
    }
}

pub struct Frequency {
    counter: Count,
    pub frequencies: Vec<f64>,
}

impl<S> Recorder<S> for Frequency
where
    S: BinaryParameterState<Move = ToggleSwap>,
{
    // Largely similar and delegated to `Count`.
    fn initialize(state: &S) -> Self {
        Self {
            counter: Count::initialize(state),
            frequencies: Vec::new(),
        }
    }

    fn record(&mut self, m: &S::Move, step: usize) {
        Recorder::<S>::record(&mut self.counter, m, step);
    }

    fn finalize(&mut self, final_step: usize) {
        Recorder::<S>::finalize(&mut self.counter, final_step);

        let total = final_step as f64;
        self.frequencies = self
            .counter
            .counts
            .iter()
            .map(|&c| c as f64 / total)
            .collect();
    }
}

impl Measure for Frequency {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.frequencies.iter().map(|&f| f as f64)
    }

    fn diagnostics(&self) -> impl Iterator<Item = Option<String>> {
        // Delegate to the counter for swap/toggle info
        self.counter.swaps.iter().map(|&x| Some(x.to_string()))
    }

    fn get_score(&self, i: usize) -> f64 {
        self.frequencies[i]
    }

    fn get_diagnostics(&self, i: usize) -> Option<f64> {
        Some(self.counter.swaps[i] as f64)
    }
}
