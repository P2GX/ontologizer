use crate::bayesian::state::{
    Increment, MgsaMove, MgsaState, ParameterState, State, TermState, ToggleSwap,
};
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

// --- Term Proposer ---

pub struct TermToggleProposer;
impl TermToggleProposer {
    pub fn new() -> Self {
        Self
    }
}

impl Proposer<TermState> for TermToggleProposer {
    /// Draw a random number representing the index to all possible toggles and swaps and propose
    /// the corresponding toggle or swap.
    fn propose<R: Rng>(&self, state: &TermState, rng: &mut R) -> <TermState as State>::Move {
        // Every possible state transition is equally likely.
        let n = state.n_terms();

        let x = rng.random_range(0..n);

        ToggleSwap::Toggle(x)
    }

    fn log_proposal(&self, state: &TermState) -> f64 {
        let n = state.n_terms();

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

impl Proposer<MgsaState> for TermToggleProposer {
    fn propose<R: Rng>(&self, state: &MgsaState, rng: &mut R) -> MgsaMove {
        let m = <Self as Proposer<TermState>>::propose(self, &state.terms, rng);
        MgsaMove::Term(m)
    }

    fn log_proposal(&self, state: &MgsaState) -> f64 {
        <Self as Proposer<TermState>>::log_proposal(self, &state.terms)
    }

    fn log_proposal_ratio(&self, state: &MgsaState, m: &MgsaMove) -> Option<f64> {
        match m {
            MgsaMove::Term(tm) => {
                <Self as Proposer<TermState>>::log_proposal_ratio(self, &state.terms, tm)
            }
            MgsaMove::Parameter(_) => None,
        }
    }
}

pub struct TermToggleSwapProposer;

impl TermToggleSwapProposer {
    pub fn new() -> Self {
        Self
    }
}

impl Proposer<TermState> for TermToggleSwapProposer {
    /// Draw a random number representing the index to all possible toggles and swaps and propose
    /// the corresponding toggle or swap.
    fn propose<R: Rng>(&self, state: &TermState, rng: &mut R) -> <TermState as State>::Move {
        // Every possible state transition is equally likely.
        let n = state.n_terms();
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
        let n = state.n_terms();
        let na = state.n_active();
        let ni = state.n_inactive();

        -((n + na * ni) as f64).ln() // It is log(1/n)
    }

    ///
    fn log_proposal_ratio(&self, state: &TermState, m: &<TermState as State>::Move) -> Option<f64> {
        match *m {
            ToggleSwap::Toggle(i) => {
                let n_current = (state.n_terms() + state.n_active() * state.n_inactive()) as f64;

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

impl Proposer<MgsaState> for TermToggleSwapProposer {
    fn propose<R: Rng>(&self, state: &MgsaState, rng: &mut R) -> <MgsaState as State>::Move {
        let m = <Self as Proposer<TermState>>::propose(self, &state.terms, rng);
        MgsaMove::Term(m)
    }

    fn log_proposal(&self, state: &MgsaState) -> f64 {
        <Self as Proposer<TermState>>::log_proposal(self, &state.terms)
    }

    fn log_proposal_ratio(&self, state: &MgsaState, m: &<MgsaState as State>::Move) -> Option<f64> {
        match m {
            MgsaMove::Term(tm) => {
                <Self as Proposer<TermState>>::log_proposal_ratio(self, &state.terms, tm)
            }
            MgsaMove::Parameter(_) => None,
        }
    }
}
// --- Parameter Proposer ---
pub struct ParameterGaussProposer {
    gaussian: Normal<f64>,
}

impl ParameterGaussProposer {
    pub fn new(sigma: f64) -> Self {
        Self {
            gaussian: Normal::new(0.0, sigma).expect("Invalid sigma for GaussianProposer"),
        }
    }
}

impl Proposer<ParameterState> for ParameterGaussProposer {
    fn propose<R: Rng>(
        &self,
        state: &ParameterState,
        rng: &mut R,
    ) -> <ParameterState as State>::Move {
        // sample on logit space ln(p / (1-p)) because [0, 1] becomes [-oo, oo]
        let index = rng.random_range(0..state.n_params());
        let x = state.get(index);
        let logit = (x / (1. - x)).ln();

        let noise = self.gaussian.sample(rng);
        let logit_new = logit + noise;

        let x_new = 1.0 / (1.0 + (-logit_new).exp());

        let delta = x_new - x;
        Increment { index, delta }
    }

    fn log_proposal(&self, state: &ParameterState) -> f64 {
        let p = state.p();
        let alpha = state.alpha();
        let beta = state.beta();

        -(p * (1. - p)).ln() - (alpha * (1. - alpha)).ln() - (beta * (1. - beta)).ln()
    }

    fn log_proposal_ratio(
        &self,
        state: &ParameterState,
        m: &<ParameterState as State>::Move,
    ) -> Option<f64> {
        let index = m.index;
        let p = state.get(index);
        let p_new = p + m.delta;

        let f1 = (p * (1.0 - p)).ln();
        let f2 = (p_new * (1.0 - p_new)).ln();

        // ln(p(1-p)) - ln(p'(1-p'))
        Some(f1 - f2)
    }
}

impl Proposer<MgsaState> for ParameterGaussProposer {
    fn propose<R: Rng>(&self, state: &MgsaState, rng: &mut R) -> <MgsaState as State>::Move {
        let m = <Self as Proposer<ParameterState>>::propose(self, &state.params, rng);
        MgsaMove::Parameter(m)
    }

    fn log_proposal(&self, state: &MgsaState) -> f64 {
        <Self as Proposer<ParameterState>>::log_proposal(self, &state.params)
    }

    fn log_proposal_ratio(&self, state: &MgsaState, m: &<MgsaState as State>::Move) -> Option<f64> {
        match m {
            MgsaMove::Term(_) => None,
            MgsaMove::Parameter(pm) => {
                <Self as Proposer<ParameterState>>::log_proposal_ratio(self, &state.params, pm)
            }
        }
    }
}
// --- Mixed Proposer ---

pub struct MixedProposer {
    term_proposer: TermToggleSwapProposer,
    param_proposer: ParameterGaussProposer,
    term_update_prob: f64,
}

impl MixedProposer {
    pub fn new(
        term_proposer: TermToggleSwapProposer,
        param_proposer: ParameterGaussProposer,
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

    fn log_proposal(&self, state: &MgsaState) -> f64 {
        self.term_proposer.log_proposal(&state.terms)
            + self.param_proposer.log_proposal(&state.params)
    }

    fn log_proposal_ratio(&self, state: &MgsaState, m: &MgsaMove) -> Option<f64> {
        match m {
            MgsaMove::Term(tm) => self.term_proposer.log_proposal_ratio(&state.terms, tm),
            MgsaMove::Parameter(pm) => self.param_proposer.log_proposal_ratio(&state.params, pm),
        }
    }
}
