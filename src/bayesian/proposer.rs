use crate::bayesian::state::{MgsaState, ParameterState, State, TermState};
use rand::{Rng, RngExt};
use rand_distr::{Distribution, Normal};

/// Allows to sample from the sample space and to find the corresponding transition probability.
pub trait Proposer<S: State> {
    /// Generates move x -> x'
    fn propose<R: Rng>(&self, state: &S, rng: &mut R) -> S::Move;
    /// Calculates q(x'|x)
    fn log_proposal(&self, state: &S) -> f64;
    /// Calculates q(x'|x) / q(x|x')
    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> Option<f64>;
}

// --- Moves ---

/// Explore the sample space by flipping the value of one element,
/// or by exchanging the (different) values of two elements.
#[derive(Clone, Copy, Debug)]
pub enum ToggleSwap {
    Toggle(usize),
    Swap(usize, usize),
}

/// A move specifically for continuous parameters.
#[derive(Clone, Copy, Debug)]
pub struct Increment {
    pub index: usize,
    pub delta: f64,
}

/// A move that can either be a Term change or a Parameter change.
#[derive(Clone, Debug)]
pub enum MgsaMove {
    Term(ToggleSwap),
    Parameter(Increment),
}

// --- Term Proposer ---

pub struct ToggleProposer;
impl ToggleProposer {
    pub fn new() -> Self {
        Self
    }
}

impl Proposer<TermState> for ToggleProposer {
    /// Draw a random number representing the index to all possible toggles and swaps and propose
    /// the corresponding toggle or swap.
    fn propose<R: Rng>(&self, state: &TermState, rng: &mut R) -> <TermState as State>::Move {
        // Every possible state transition is equally likely.
        let n = state.n_all();

        let x = rng.random_range(0..n);

        ToggleSwap::Toggle(x)
    }

    fn log_proposal(&self, state: &TermState) -> f64 {
        let n = state.n_all();

        -(n as f64).ln() // It is log(1/n)
    }

    ///
    fn log_proposal_ratio(
        &self,
        _state: &TermState,
        _m: &<TermState as State>::Move,
    ) -> Option<f64> {
        Some(0.0)
    }
}

pub struct ToggleSwapProposer;
impl ToggleSwapProposer {
    pub fn new() -> Self {
        Self
    }
}

impl Proposer<TermState> for ToggleSwapProposer {
    /// Draw a random number representing the index to all possible toggles and swaps and propose
    /// the corresponding toggle or swap.
    fn propose<R: Rng>(&self, state: &TermState, rng: &mut R) -> <TermState as State>::Move {
        // Every possible state transition is equally likely.
        let n = state.n_all();
        let na = state.n_active();
        let ni = state.n_inactive();

        let m = n + na * ni;
        let x = rng.random_range(0..m);

        // In m cases we flip the state of a single term.
        if x < n {
            ToggleSwap::Toggle(x)
        }
        // in m_on * m_off cases we swap the states of two terms with different states.
        else {
            // map random number x to pairs of indices a, b
            let k = x - n;
            let i = state.get_active(k % na);
            let j = state.get_inactive(k / na);
            ToggleSwap::Swap(i, j)
        }
    }

    fn log_proposal(&self, state: &TermState) -> f64 {
        let n = state.n_all();
        let na = state.n_active();
        let ni = state.n_inactive();

        -((n + na * ni) as f64).ln() // It is log(1/n)
    }

    ///
    fn log_proposal_ratio(&self, state: &TermState, m: &<TermState as State>::Move) -> Option<f64> {
        match *m {
            ToggleSwap::Toggle(i) => {
                let n_current = (state.n_all() + state.n_active() * state.n_inactive()) as f64;

                // Calculate the change in possible moves (delta).

                // Active -> Inactive: new_n = current_n + n_on - n_off - 1
                // Inactive -> Active: new_n = current_n + n_off - n_on - 1
                let diff = state.n_active() as f64 - state.n_inactive() as f64;
                // If terms[i] is true (Active->Inactive), we use 'diff-1'.
                // If terms[i] is false (Inactive->Active), we use '-diff-1'.
                let delta = if state.get(i) { diff } else { -diff } - 1.0;

                let n_proposed = n_current + delta;
                Some((n_proposed / n_current).ln())
            }
            // Swapping preserves n_on and n_off, so the state space size N stays constant.
            // ln(N / N) = ln(1) = 0
            ToggleSwap::Swap(_, _) => Some(0.0),
        }
    }
}

// --- Parameter Proposer ---
pub struct GaussianProposer {
    distribution: Normal<f64>,
}

impl GaussianProposer {
    pub fn new(sigma: f64) -> Self {
        Self {
            distribution: Normal::new(0.0, sigma).expect("Invalid sigma for GaussianProposer"),
        }
    }
}

impl Proposer<ParameterState> for GaussianProposer {
    fn propose<R: Rng>(
        &self,
        state: &ParameterState,
        rng: &mut R,
    ) -> <ParameterState as State>::Move {
        let index = rng.random_range(0..state.n_all());
        let old_value = state.get(index);

        let noise = self.distribution.sample(rng);
        let mut new_value = old_value + noise;

        // Reflection for [0, 1] bounds with loop to ensure safety
        while new_value < 0.0 || new_value > 1.0 {
            if new_value < 0.0 {
                new_value = -new_value;
            } else if new_value > 1.0 {
                new_value = 2.0 - new_value;
            }
        }

        // We store the effective delta so apply/revert is simple arithmetic
        let delta = new_value - old_value;

        Increment { index, delta }
    }

    fn log_proposal(&self, _state: &ParameterState) -> f64 {
        0.0 // Symmetric proposal
    }

    fn log_proposal_ratio(
        &self,
        _state: &ParameterState,
        _m: &<ParameterState as State>::Move,
    ) -> Option<f64> {
        Some(0.0) // Symmetric proposal (Gaussian)
    }
}

// --- Mixed Proposer ---

pub struct MixedProposer {
    term_proposer: ToggleSwapProposer,
    param_proposer: GaussianProposer,
    term_update_prob: f64,
}

impl MixedProposer {
    pub fn new(
        term_proposer: ToggleSwapProposer,
        param_proposer: GaussianProposer,
        term_update_prob: f64,
    ) -> Self {
        Self {
            term_proposer,
            param_proposer,
            term_update_prob,
        }
    }
}

impl Proposer<MgsaState> for MixedProposer {
    fn propose<R: Rng>(&self, state: &MgsaState, rng: &mut R) -> <MgsaState as State>::Move {
        if rng.random_bool(self.term_update_prob) {
            let m = self.term_proposer.propose(&state.terms, rng);
            MgsaMove::Term(m)
        } else {
            let m = self.param_proposer.propose(&state.params, rng);
            MgsaMove::Parameter(m)
        }
    }

    fn log_proposal(&self, _state: &MgsaState) -> f64 {
        0.0
    }

    fn log_proposal_ratio(&self, state: &MgsaState, m: &MgsaMove) -> Option<f64> {
        match m {
            MgsaMove::Term(tm) => self.term_proposer.log_proposal_ratio(&state.terms, tm),
            MgsaMove::Parameter(pm) => self.param_proposer.log_proposal_ratio(&state.params, pm),
        }
    }
}
