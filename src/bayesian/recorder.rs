use crate::bayesian::measure::Probability;
use crate::bayesian::proposer::ToggleSwap;
use crate::bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::bayesian::state::{BinaryParameterState, State};
pub trait Recorder<S: State> {
    type Target;
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, m: &<S as State>::Move, step: usize);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(self, final_step: usize) -> Self::Target;
}


pub struct ProbabilityRecorder {
    counts: Vec<usize>,
    swaps: Vec<usize>,
    active_since: Vec<Option<usize>>,
}

impl ProbabilityRecorder {
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

impl<S> Recorder<S> for ProbabilityRecorder
where
    S: BinaryParameterState<Move = ToggleSwap>,
{
    type Target = Probability;
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

    fn finalize(mut self, final_step: usize) -> Self::Target {
        let mut probabilities = vec![0.0; self.counts.len()];
        for (i, start_opt) in self.active_since.iter().enumerate() {
            if let Some(start) = start_opt {
                self.counts[i] += final_step - start;
            }
            probabilities[i] = self.counts[i] as f64 / final_step as f64;
        }

        Probability::new(probabilities, self.swaps)
    }
}

