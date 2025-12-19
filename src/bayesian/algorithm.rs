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

pub struct MetropolisHasting<M: Model, P> {
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

    fn get_log_proposal_ratio(
        &mut self,
        state: &mut M::State,
        m: &<M::State as State>::Move,
    ) -> f64 {
        match self.proposer.log_proposal_ratio(state, &m) {
            Some(log_p_ratio) => log_p_ratio,
            None => {
                let log_p1 = self.proposer.log_proposal(state);
                state.apply(&m);
                let log_p2 = self.proposer.log_proposal(state);
                state.revert(&m);
                log_p2 - log_p1
            }
        }
    }

    fn get_log_prior_ratio(&mut self, state: &mut M::State, m: &<M::State as State>::Move) -> f64 {
        match self.model.log_prior_ratio(state, &m) {
            Some(log_p_ratio) => log_p_ratio,
            None => {
                let log_p1 = self.model.log_prior(state);
                state.apply(&m);
                let log_p2 = self.model.log_prior(state);
                state.revert(&m);
                log_p2 - log_p1
            }
        }
    }

    fn get_log_likelihood_ratio(
        &mut self,
        state: &mut M::State,
        cache: &M::Cache,
        m: &<M::State as State>::Move,
    ) -> f64 {
        match self.model.log_likelihood_ratio(cache, m) {
            Some(log_l_ratio) => log_l_ratio,
            None => {
                let log_l1 = self.model.log_likelihood(cache);
                state.apply(&m);
                let log_l2 = self.model.log_likelihood(cache);
                state.revert(&m);
                log_l2 - log_l1
            }
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

        let mut cache = self.model.create_cache(state);

        let mut rng = rand::rng();

        let n_max = self.burn_in + self.iterations;
        for i in 0..n_max {
            let m = self.proposer.propose(state, &mut rng);

            let log_q_ratio = self.get_log_proposal_ratio(state, &m);

            let log_p_ratio = self.get_log_prior_ratio(state, &m);

            let log_l_ratio = self.get_log_likelihood_ratio(state, &cache, &m);

            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let x: f64 = rng.random_range(0.0..1.0);
            if log_accept >= 0.0 || x.ln() < log_accept {
                state.apply(&m);

                self.model.update_cache(&mut cache, state, &m);

                if i > self.burn_in {
                    result.record(&m, i);
                }
            }
        }
        result.finalize(n_max);
        result
    }
}
