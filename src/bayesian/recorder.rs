use crate::bayesian::measure::{Mean, Probability};
use crate::bayesian::state::{MgsaMove, MgsaState, State, ToggleSwap};

pub trait Recorder<S: State> {
    type Target;
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, m: &<S as State>::Move, step: usize);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(self, final_step: usize) -> Self::Target;
}

pub struct TermRecorder {
    counts: Vec<usize>,
    swaps: Vec<usize>,
    active_since: Vec<Option<usize>>,
}

impl TermRecorder {
    fn update_counts(&mut self, idx: usize, step: usize) {
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

impl Recorder<MgsaState> for TermRecorder {
    type Target = Vec<Probability>;
    fn initialize(state: &MgsaState) -> Self {
        let n = state.terms.n_terms();
        let mut active_since = vec![None; n];

        for i in 0..n {
            if state.terms.get(i) {
                active_since[i] = Some(0);
            }
        }

        Self {
            counts: vec![0; n],
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
            MgsaMove::Parameter(_) => {}
        }
    }

    fn finalize(mut self, final_step: usize) -> Self::Target {
        let mut probabilities = Vec::new();
        for (i, start_opt) in self.active_since.iter().enumerate() {
            if let Some(start) = start_opt {
                self.counts[i] += final_step - start;
            }
            probabilities.push(Probability::new(
                self.counts[i] as f64 / final_step as f64,
                self.swaps[i],
            ))
        }

        probabilities
    }
}

// --- Parameter Recorder ---
pub struct ParameterRecorder {
    current_values: [f64; 3],
    sums: [f64; 3],
    changes: [usize; 3],
    last_change_step: usize,
}

impl ParameterRecorder {
    fn update_sums(&mut self, step: usize) {
        let duration = (step - self.last_change_step) as f64;
        for i in 0..3 {
            self.sums[i] += self.current_values[i] * duration;
        }
        self.last_change_step = step;
    }
}

impl Recorder<MgsaState> for ParameterRecorder {
    type Target = Vec<Mean>;

    fn initialize(state: &MgsaState) -> Self {
        Self {
            current_values: [state.params.p(), state.params.alpha(), state.params.beta()],
            sums: [0.0; 3],
            changes: [0; 3],
            last_change_step: 0,
        }
    }

    fn record(&mut self, m: &<MgsaState as State>::Move, step: usize) {
        match m {
            MgsaMove::Term(_) => {}
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

// --- Composite Recorder ---
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
