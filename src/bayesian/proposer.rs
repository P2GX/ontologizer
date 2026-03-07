use crate::bayesian::state::{
    Increment, MgsaMove, MgsaState, ParameterState, State, TermState, ToggleSwap,
};
use rand::{Rng, RngExt};
use rand_distr::{Distribution, Normal};

/// A generic interface for generating transitions (moves) in the state space.
///
/// It handles the proposal distribution $q(x'|x)$ required for the Metropolis-Hastings algorithm.
pub trait Proposer<S: State> {
    /// Generates move x -> x'
    fn propose<R: Rng>(&self, state: &S, rng: &mut R) -> S::Move;

    /// Calculates log q(x'|x), up to an additive constant.
    ///
    /// Constant factors that are identical for the forward and reverse proposals — such as
    /// normalisation constants, symmetric density terms, or terms depending only on unchanged
    /// parts of the state — may be omitted, because this value is only used inside
    /// `log_proposal_ratio` where such terms cancel.
    fn log_proposal(&self, state: &S) -> f64;

    /// Calculates log( q(x'|x) / q(x|x') )
    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> Option<f64>;
}

// ==========================================
// TERM PROPOSERS
// ==========================================

/// A simple proposer that only toggles terms active/inactive.
///
/// It chooses a term uniformly at random and flips its state.
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
        // Uniformly select a term inde
        let n_terms = state.n_terms();
        let idx = rng.random_range(0..n_terms);
        ToggleSwap::Toggle(idx)
    }

    fn log_proposal(&self, state: &TermState) -> f64 {
        // log q(x'|x) = -ln(N). Exact — no constant factors omitted.
        let n_terms = state.n_terms();
        -(n_terms as f64).ln()
    }

    ///
    fn log_proposal_ratio(
        &self,
        _state: &TermState,
        _m: &<TermState as State>::Move,
    ) -> Option<f64> {
        // The distribution is symmetric: q(x'|x) = q(x|x') = 1/N
        // log(1) = 0
        Some(0.0)
    }
}

// Implement for MgsaState wrapper
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

/// A proposer that mixes simple Toggles with Swaps.
///
/// A "Swap" exchanges an active term with an inactive one, preserving the total count.
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
        let n_terms = state.n_terms();
        let n_active = state.n_active();
        let n_inactive = state.n_inactive();

        let m = n_terms + n_active * n_inactive;
        let choice = rng.random_range(0..m);

        if choice < n_terms {
            // Case 1: Toggle
            ToggleSwap::Toggle(choice)
        } else {
            // Case 2: Swap
            // Map the random choice back to a specific (Active, Inactive) pair
            let k = choice - n_terms;
            let active_idx = state.get_active(k % n_active);
            let inactive_idx = state.get_inactive(k / n_active);
            ToggleSwap::Swap(active_idx, inactive_idx)
        }
    }

    fn log_proposal(&self, state: &TermState) -> f64 {
        // log q(x'|x) = -ln(M). Exact — no constant factors omitted.
        let n_terms = state.n_terms();
        let n_active = state.n_active();
        let n_inactive = state.n_inactive();

        let m = n_terms + n_active * n_inactive;
        -(m as f64).ln()
    }

    ///
    fn log_proposal_ratio(&self, state: &TermState, m: &<TermState as State>::Move) -> Option<f64> {
        match *m {
            ToggleSwap::Toggle(i) => {
                let n_terms = state.n_terms();
                let n_active = state.n_active();
                let n_inactive = state.n_inactive();
                let current_space_size = (n_terms + n_active * n_inactive) as f64;

                // Calculate the new space size after the toggle
                // Active -> Inactive: n_active decreases by 1, n_inactive increases by 1
                // Inactive -> Active: n_active increases by 1, n_inactive decreases by 1
                let diff = state.n_active() as f64 - state.n_inactive() as f64;

                // If terms[i] is true (Active->Inactive), new product is (na-1)(ni+1) = na*ni + na - ni - 1
                // Delta = na - ni - 1

                // If terms[i] is false (Inactive->Active), new product is (na+1)(ni-1) = na*ni - na + ni - 1
                // Delta = -na + ni - 1 = -(na - ni) - 1

                let delta = if state.get(i) { diff } else { -diff } - 1.0;
                let proposed_space_size = current_space_size + delta;

                // log( q(x'|x) / q(x|x')
                // Ratio = (1/current) / (1/proposed) = proposed / current.
                // Log = ln(proposed) - ln(current).
                Some((proposed_space_size / current_space_size).ln())
            }
            // Swapping preserves na and ni, so the state space size N stays constant.
            // ln(N / N) = ln(1) = 0
            ToggleSwap::Swap(_, _) => Some(0.0),
        }
    }
}

// Implement for MgsaState wrapper
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

// ==========================================
// PARAMETER PROPOSERS
// ==========================================

/// A Random Walk Metropolis proposer for continuous parameters (p, alpha, beta).
///
/// It proposes moves in the logit space to respect the [0, 1] bounds naturally,
/// adding Gaussian noise to the transformed variable.
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
        // Pick a parameter index uniformly
        let idx = rng.random_range(0..state.n_params());
        let parameter = state.get(idx);

        // Transform to logit space: y = ln(x / (1-x))
        let logit = (parameter / (1. - parameter)).ln();

        // Add Gaussian noise: N(0, sigma)
        let noise = self.gaussian.sample(rng);
        let logit_new = logit + noise;

        // Transform back: x' = sigmoid(y')
        let x_new = 1.0 / (1.0 + (-logit_new).exp());

        let delta = x_new - parameter;
        Increment { index: idx, delta }
    }

    fn log_proposal(&self, state: &ParameterState) -> f64 {
        // log q(x'|x) up to an additive constant. The omitted terms are:
        //   - ln(1/3): uniform probability of selecting a parameter (same in forward and reverse)
        //   - log φ(logit(x') - logit(x)): Gaussian density, cancels due to symmetry φ(a) = φ(-a)
        //   - Jacobian terms for the two unchanged parameters (identical in both states)
        // Only the Jacobian of the changed parameter survives in the ratio, so the sum over
        // all three Jacobians here is sufficient when used exclusively in log_proposal_ratio.
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
        // q(x'|x) for a logit-normal random walk includes the Jacobian of the transformation.
        // Ratio = (x(1-x)) / (x'(1-x'))
        // Log Ratio = ln(x(1-x)) - ln(x'(1-x'))
        let idx = m.index;
        let parameter = state.get(idx);
        let parameter_new = parameter + m.delta;

        let log_jac_current = (parameter * (1.0 - parameter)).ln();
        let log_jac_proposed = (parameter_new * (1.0 - parameter_new)).ln();

        // ln(p(1-p)) - ln(p'(1-p'))
        Some(log_jac_current - log_jac_proposed)
    }
}

// Implement for MgsaState wrapper
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

// ==========================================
// MIXED PROPOSER
// ==========================================

/// A proposer that selects between Term moves and Parameter moves based on a probability.
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
        // log q(x'|x) up to an additive constant. The mixing weight (term_update_prob) is omitted
        // because it is identical for forward and reverse proposals and cancels in the ratio.
        // For a term move the parameter Jacobian terms are the same in both states and cancel;
        // for a parameter move the term space-size term is the same in both states and cancels.
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
