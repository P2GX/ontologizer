use crate::better_bayesian::observation::Observation;
use crate::better_bayesian::state::{CountableState, State};

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model<S: CountableState, O: Observation>{

    // Log probability P(S) to find a State configuration
    fn log_prior(&self, state: &S) -> f64;

    // Log probability P(O | S) to find an Observable configuration given a State configuration.
    fn log_likelihood(&self, state : &S, observation : &O) -> f64;

    // TODO: I am not sure if we can calculate this. But it would be very useful.
    fn log_likelihood_ratio(&self, state : &S, observation : &O, m : S::Move) -> f64;

}

struct OrModel{
    latent_state : Vec<bool>,
    p : f64, // probability of term being active
    alpha : f64, // probability of a gene being incorrectly active.
    beta : f64 // The probability of a gene being incorrectly inactive.
}

impl<S : CountableState, O : Observation> Model<S, O> for OrModel {

    // Log probability P(T) to find a Term configuration
    fn log_prior(&self, state : &S) -> f64{
        let m0 = state.n_active() as f64;
        let m1 =  state.n_inactive() as f64;

        m0 * self.p.ln() + m1 * (1. - self.p).ln()
    }

    // Log probability P(O | T) to find an observed Gene configuration given a Terms configuration.
    fn log_likelihood(&self, state: &S, observation: &O) -> f64 {
        todo!()
    }

    fn log_likelihood_ratio(&self, state: &S, observation: &O, m: S::Move) -> f64 {
        todo!()
    }
}