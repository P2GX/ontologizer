use rand::Rng;
use crate::better_bayesian::state::State;

/// Represents the (Marginal) Probability P(x).
pub trait Probability {
    type Event: State;

    fn probability(&self, x: &Self::Event) -> f64 {
        self.log_probability(x).exp()
    }

    fn log_probability(&self, x: &Self::Event) -> f64;
}

/// Represents the Conditional or Transition Probability P(x_{t+1} | x_t).
pub trait TransitionProbability {
    type Event: State;

    fn conditional_prob(&self, target: &Self::Event, source: &Self::Event) -> f64 {
        self.log_conditional_prob(target, source).exp()
    }

    fn log_conditional_prob(&self, target: &Self::Event, source: &Self::Event) -> f64;
}