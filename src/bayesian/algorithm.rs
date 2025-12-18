// MCMC
use crate::bayesian::model::Model;
use crate::bayesian::proposer::Proposer;
use crate::bayesian::recorder::Recorder;
use crate::bayesian::state::State;
use rand::Rng;

pub trait Algorithm<M>
where
    M: Model,
{
    fn sample<R: Recorder<M::State>>(&mut self, state: &mut M::State) -> R;
}

pub struct MetropolisHasting<M, P> {
    model: M,
    proposer: P,
    iterations: usize,
    burn_in: usize,
}

impl<M, P> MetropolisHasting<M, P>
where
    M: Model,
    P: Proposer<M::State>,
{
    pub fn new(model: M, proposer: P, iterations: usize, burn_in: usize) -> Self {
        Self {
            model,
            proposer,
            iterations,
            burn_in,
        }
    }
}

impl<M, P> Algorithm<M> for MetropolisHasting<M, P>
where
    M: Model,
    P: Proposer<M::State>,
{
    fn sample<R>(&mut self, state: &mut M::State) -> R
    where
        R: Recorder<M::State>,
    {
        let mut result = R::initialize(&state);

        let mut rng = rand::rng();

        for i in 0..(self.burn_in + self.iterations) {
            if i % 1000 == 0 {
                println!("{i:?}")
            }
            let m = self.proposer.propose(state, &mut rng);
            let log_q_ratio = self.proposer.log_proposal_ratio(state, &m);

            let log_p_ratio = match self.model.log_prior_ratio(state, &m) {
                Some(log_p_ratio) => log_p_ratio,
                None => {
                    let log_q1 = self.model.log_prior(state);
                    state.apply(&m);
                    let log_q2 = self.model.log_prior(state);
                    state.revert(&m);
                    log_q2 - log_q1
                }
            };

            let log_l_ratio = match self.model.log_likelihood_ratio(state, &m) {
                Some(log_p_ratio) => log_p_ratio,
                None => {
                    let log_l1 = self.model.log_likelihood(state);
                    state.apply(&m);
                    let log_l2 = self.model.log_likelihood(state);
                    state.revert(&m);
                    log_l2 - log_l1
                }
            };

            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let x: f64 = rng.random_range(0.0..1.0);
            if log_accept >= 0.0 || x.ln() < log_accept {
                state.apply(&m);
            }

            if i > self.burn_in {
                result.record(&state);
            }
        }
        result.finalize();
        result
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_algorithm() {
        todo!()
    }
}
