use crate::better_bayesian::observation::{Genes, Observation};
use crate::better_bayesian::proposer::ToggleSwap;
use crate::better_bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::better_bayesian::state::{CountableState, State};
use crate::core::AnnotationIndex;
use std::sync::Arc;

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model {
    type State: State;
    type Observation: Observation;
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
    fn log_likelihood(&self, state: &Self::State, observation: &Self::Observation) -> f64;

    // Fast implementation for log likelihood ratio P(O|S2)/P(O|S1). May not be available for the specific move.
    fn log_likelihood_ratio(
        &self,
        state: &Self::State,
        observation: &Self::Observation,
        m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }
}

struct OrModel<S> {
    annotations: Arc<AnnotationIndex>,
    p: f64,     // probability of term being active
    alpha: f64, // probability of a gene being incorrectly active.
    beta: f64,  // The probability of a gene being incorrectly inactive.
    _marker: std::marker::PhantomData<S>,
}

impl<S> Model for OrModel<S>
where
    S: CountableState<Move = ToggleSwap>,
{
    type State = S;
    type Observation = Genes;

    // Log probability P(T) to find a Term configuration
    fn log_prior(&self, state: &S) -> f64 {
        let m0 = state.n_active() as f64;
        let m1 = state.n_inactive() as f64;

        m0 * self.p.ln() + m1 * (1. - self.p).ln()
    }

    fn log_prior_ratio(&self, state: &Self::State, m: &ToggleSwap) -> Option<f64> {
        let ratio = match *m {
            Toggle(i) => {
                if state.get(i) {
                    ((1.0 - self.p) / (self.p)).ln()
                } else {
                    (self.p / (1.0 - self.p)).ln()
                }
            }
            Swap(_, _) => 0.0,
        };
        Some(ratio)
    }

    // Log probability P(O | T) to find an observed Gene configuration given a Terms configuration.
    fn log_likelihood(&self, state: &S, observation: &Genes) -> f64 {
        todo!()
    }

    fn log_likelihood_ratio(&self, state: &S, observation: &Genes, m: &ToggleSwap) -> Option<f64> {
        todo!()
    }
}
