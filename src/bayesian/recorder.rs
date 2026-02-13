use crate::bayesian::measure::{Mean, Probability};
use crate::bayesian::state::{MgsaMove, MgsaState, State, ToggleSwap};

/// A generic interface for recording the history of a Markov Chain.
pub trait Recorder<S: State> {
    type Target;

    /// Initializes the recorder based on the initial state.
    fn initialize(state: &S) -> Self;

    /// Records a single step in the chain.
    fn record(&mut self, m: &<S as State>::Move, step: usize);

    /// Finalizes the recording process and produces the result.
    fn finalize(self, final_step: usize) -> Self::Target;
}

// ==========================================
// TERM RECORDER
// ==========================================

/// Records the posterior probabilities of terms being active.
pub struct TermRecorder {
    active_count: Vec<usize>,
    swaps: Vec<usize>,
    active_since: Vec<Option<usize>>,
}

impl TermRecorder {
    /// Updates the internal timers for a specific term index.
    fn update_counts(&mut self, idx: usize, step: usize) {
        self.swaps[idx] += 1;
        match self.active_since[idx] {
            Some(start) => {
                // Was ON, turning OFF. Record the duration.
                self.active_count[idx] += step - start;
                self.active_since[idx] = None;
            }
            None => {
                // Was OFF, turning ON. Start the timer.
                self.active_since[idx] = Some(step);
            }
        }
    }
}

impl Recorder<MgsaState> for TermRecorder {
    type Target = Vec<Probability>;
    fn initialize(state: &MgsaState) -> Self {
        let n = state.terms.n_terms();
        let mut active_since = vec![None; n];

        // Initialize timers for terms that start active
        for i in 0..n {
            if state.terms.get(i) {
                active_since[i] = Some(0);
            }
        }

        Self {
            active_count: vec![0; n],
            swaps: vec![0; n],
            active_since,
        }
    }

    fn record(&mut self, m: &<MgsaState as State>::Move, step: usize) {
        match m {
            MgsaMove::Term(tm) => match tm {
                ToggleSwap::Toggle(i) => self.update_counts(*i, step),
                ToggleSwap::Swap(i, j) => {
                    self.update_counts(*i, step);
                    self.update_counts(*j, step);
                }
            },
            MgsaMove::Parameter(_) => {
                // Parameter moves do not affect term states.
            }
        }
    }

    fn finalize(mut self, final_step: usize) -> Self::Target {
        let mut result = Vec::new();
        for (i, start_opt) in self.active_since.iter().enumerate() {
            if let Some(start) = start_opt {
                self.active_count[i] += final_step - start;
            }
            result.push(Probability::new(
                self.active_count[i] as f64 / final_step as f64,
                self.swaps[i],
            ))
        }

        result
    }
}

// ==========================================
// PARAMETER RECORDER
// ==========================================

/// Records the posterior means of continuous parameters (p, alpha, beta).

pub struct ParameterRecorder {
    current_values: [f64; 3],
    sums: [f64; 3],
    changes: [usize; 3],
    last_update_step: usize,
}

impl ParameterRecorder {
    /// Adds the contribution of the current values to the sums for the duration since the last update.
    fn update_sums(&mut self, step: usize) {
        let duration = (step - self.last_update_step) as f64;
        for i in 0..3 {
            self.sums[i] += self.current_values[i] * duration;
        }
        self.last_update_step = step;
    }
}

impl Recorder<MgsaState> for ParameterRecorder {
    type Target = Vec<Mean>;

    fn initialize(state: &MgsaState) -> Self {
        Self {
            current_values: [state.params.p(), state.params.alpha(), state.params.beta()],
            sums: [0.0; 3],
            changes: [0; 3],
            last_update_step: 0,
        }
    }

    fn record(&mut self, m: &<MgsaState as State>::Move, step: usize) {
        match m {
            MgsaMove::Term(_) => {
                // Term moves do not change parameters.
            }
            MgsaMove::Parameter(pm) => {
                self.update_sums(step);
                if pm.index < 3 {
                    self.current_values[pm.index] += pm.delta;
                    self.changes[pm.index] += 1;
                }
            }
        }
    }

    fn finalize(mut self, final_step: usize) -> Self::Target {
        self.update_sums(final_step);
        let mut means = Vec::with_capacity(3);
        for i in 0..3 {
            means.push(Mean::new(self.sums[i] / final_step as f64, self.changes[i]));
        }
        means
    }
}

// ==========================================
// COMPOSITE RECORDER
// ==========================================

pub struct MgsaRecorder {
    term_recorder: TermRecorder,
    param_recorder: ParameterRecorder,
}

impl Recorder<MgsaState> for MgsaRecorder {
    type Target = (Vec<Probability>, Vec<Mean>);

    fn initialize(state: &MgsaState) -> Self {
        Self {
            term_recorder: TermRecorder::initialize(&state),
            param_recorder: ParameterRecorder::initialize(&state),
        }
    }

    fn record(&mut self, m: &MgsaMove, step: usize) {
        self.term_recorder.record(m, step);
        self.param_recorder.record(m, step);
    }

    fn finalize(self, final_step: usize) -> Self::Target {
        (
            self.term_recorder.finalize(final_step),
            self.param_recorder.finalize(final_step),
        )
    }
}
