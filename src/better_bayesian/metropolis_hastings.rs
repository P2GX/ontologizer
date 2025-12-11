// MCMC
use rand::Rng;
use crate::better_bayesian::probability::{TransitionProbability, Probability};
use crate::better_bayesian::state::State;

pub struct MetropolisHasting<'s, P, Q, S, R>
where
    S: State,
    P: Probability<Event = S>,
    Q: TransitionProbability<Event = S>,
    R: Rng
{
    target_prob: P,
    transition_prob: Q,
    current_state: &'s S,
    rng : R
}

impl<'s, P, Q, S, R> MetropolisHasting<'s, P, Q, S, R>
where
    S: State,
    P: Probability<Event = S>,
    Q: TransitionProbability<Event = S>,
    R: Rng
{
    fn new(target_prob : P, transition_prob : Q, current_state : &'s S, rng : R) -> Self{
        Self {
            target_prob,
            transition_prob,
            current_state,
            rng
        }
    }

    fn accept_log_ratio(&self, proposed_state : &S) -> f64 {
        let log_q_forward = self.transition_prob.log_conditional_prob(&self.current_state, &proposed_state);
        let log_q_backwards = self.transition_prob.log_conditional_prob(&proposed_state, &self.current_state);

        let log_p_current = self.target_prob.log_probability(&self.current_state);
        let log_p_proposed = self.target_prob.log_probability(proposed_state);

        log_p_proposed + log_q_backwards - log_p_current - log_q_forward
    }
}


trait Operator<'s, S, R>
where
    R: Rng,
    S: State,
{
    fn advance(&mut self, state : &mut S);
}


impl<'s, P, Q, S, R> Operator<'s, S, R> for MetropolisHasting<'s, P, Q, S, R>
where
    S: State,
    P: Probability<Event = S>,
    Q: TransitionProbability<Event = S>,
    R: Rng
{
    fn advance(&mut self, state : &mut S) {

        let m = state.draw_move(&mut self.rng);
        let proposed_state = state.apply(&m);

        let log_alpha = self.accept_log_ratio(&state);

        // log_alpha >= 0 -> alpha >= 1 -> always accept.
        let accept = if log_alpha >= 0.0 {
            true
        } else {
            // a uniformly distributed value between 0 and 1:
            let x: f64 = self.rng.random();
            let log_x = x.ln();
            log_x < log_alpha
        };

        // if accept {
        //     self.current_state = state
        // }

    }
}

