use crate::bayesian::observation::{Observation};
use crate::bayesian::proposer::ToggleSwap;
use crate::bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::bayesian::state::{CountableState, State};
use crate::core::AnnotationIndex;
use rand::Rng;

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model{
    type State: State;
    type Observation: Observation;

    // Initialize the State
    fn initialize_state<R: Rng>(&self, rng: &mut R) -> &Self::State;

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

pub struct OrModel<'a, S, O> {
    annotations: &'a AnnotationIndex, // provides terms_to_genes and genes_to_terms.
    p: f64,     // probability of term being active
    alpha: f64, // probability of a gene being incorrectly observed active.
    beta: f64,  // probability of a gene being incorrectly observed inactive.
    _marker: std::marker::PhantomData<(S, O)>,
}

impl<'a, S, O> OrModel<'a, S, O>
{
    pub fn new(annotations: &'a AnnotationIndex, p: f64, alpha: f64, beta: f64) -> OrModel<'a, S, O> {
        Self {
            annotations,
            p,
            alpha,
            beta,
            // Hidden implementation detail
            _marker: std::marker::PhantomData,
        }
    }
}

impl<'a, S, O> Model for OrModel<'a, S, O>
where
    S: CountableState<Move = ToggleSwap>,
    O: Observation
{
    type State = S;
    type Observation = O;

    fn initialize_state<R: Rng>(&self, rng: &mut R) -> &Self::State {
        todo!()


    }

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
    fn log_likelihood(&self, state: &S, observation: &Self::Observation) -> f64 {
        todo!()
    }

    fn log_likelihood_ratio(&self, state: &S, observation: &Self::Observation, m: &ToggleSwap) -> Option<f64> {
        todo!()
    }
}


#[cfg(test)]
mod test {

    fn test_model(){

    }
}
