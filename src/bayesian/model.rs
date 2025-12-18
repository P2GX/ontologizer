use crate::bayesian::proposer::ToggleSwap;
use crate::bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::bayesian::state::{AlignedBinaryLatentState, BinaryParameterState, State};
use rand::Rng;

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model {
    type State: State;

    // Initialize the State
    fn sample_prior<R: Rng>(&self, rng: &mut R, n: usize) -> Vec<<Self::State as State>::Value>;

    // Log probability P(S) to find a State configuration
    fn log_prior(&self, state: &Self::State) -> f64;

    // Fast implementation for log probability ratio P(S2)/P(S1). May not be available for the specific move.
    fn log_prior_ratio(
        &self,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }

    // Log probability P(O | S) to find an Observable configuration given a State configuration.
    fn log_likelihood(&self, state: &Self::State) -> f64;

    // Fast implementation for log likelihood ratio P(O|S2)/P(O|S1). May not be available for the specific move.
    fn log_likelihood_ratio(
        &self,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }
}

pub struct OrModel<S> {
    terms_to_genes: Vec<Vec<usize>>, // provides terms_to_genes and genes_to_terms.
    observe_genes: Vec<bool>,
    p: f64,     // probability of term being active
    alpha: f64, // probability of a gene being incorrectly observed active.
    beta: f64,  // probability of a gene being incorrectly observed inactive.
    _marker: std::marker::PhantomData<S>,
}

impl<S> OrModel<S> {
    pub fn new(
        terms_to_genes: Vec<Vec<usize>>,
        observe_genes: Vec<bool>,
        p: f64,
        alpha: f64,
        beta: f64,
    ) -> OrModel<S> {
        Self {
            terms_to_genes,
            observe_genes,
            p,
            alpha,
            beta,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<S> Model for OrModel<S>
where
    S: BinaryParameterState<Move = ToggleSwap> + AlignedBinaryLatentState,
{
    type State = S;

    fn sample_prior<R: Rng>(&self, rng: &mut R, n_terms: usize) -> Vec<bool> {
        (0..n_terms).map(|_| rng.random_bool(self.p)).collect()
    }

    // Log probability P(T) to find a Term configuration
    fn log_prior(&self, state: &S) -> f64 {
        let m0 = state.n_active() as f64;
        let m1 = state.n_inactive() as f64;

        m0 * self.p.ln() + m1 * (1. - self.p).ln()
    }

    fn log_prior_ratio(&self, state: &Self::State, m: &ToggleSwap) -> Option<f64> {
        let log_ratio = match *m {
            Toggle(i) => {
                if state.get(i) {
                    // from active to inactive
                    ((1.0 - self.p) / (self.p)).ln()
                } else {
                    // from inactive to active
                    (self.p / (1.0 - self.p)).ln()
                }
            }
            Swap(_, _) => 0.0,
        };
        Some(log_ratio)
    }

    // Log probability P(O | T) to find an observed Gene configuration given a Terms configuration.
    fn log_likelihood(&self, state: &S) -> f64 {
        let m00 = state.n_true_positive() as f64;
        let m01 = state.n_false_positive() as f64;
        let m11 = state.n_true_negative() as f64;
        let m10 = state.n_false_negative() as f64;

        m00 * (1. - self.alpha).ln()
            + m01 * self.alpha.ln()
            + m11 * (1. - self.beta).ln()
            + m10 * self.beta.ln()
    }

    fn log_likelihood_ratio(&self, state: &S, m: &ToggleSwap) -> Option<f64> {
        None
    }
}
