// MCMC
use crate::better_bayesian::model::Model;
use crate::better_bayesian::proposer::Proposer;
use crate::better_bayesian::state::{CountableState, State};
use rand::Rng;

pub trait Recorder<S: State> {
    /// Initialize the recorder (e.g., allocate vectors based on state size)
    fn initialize(state: &S) -> Self;

    /// Record a single step
    fn record(&mut self, state: &S);

    // Optional: A method to finalize computation (e.g., divide sum by count)
    fn finalize(&mut self);
}

struct Count {
    counts: Vec<u32>,
    n: u32,
}
impl<S> Recorder<S> for Count
where
    S: CountableState,
{
    fn initialize(state: &S) -> Self {
        Self {
            counts: vec![0; state.n_all()],
            n: 0,
        }
    }
    fn record(&mut self, state: &S) {
        for i in 0..state.n_all() {
            if state.get(i) {
                self.counts[i] += 1;
            }
        }
        self.n += 1;
    }

    fn finalize(&mut self) {}
}

struct Frequency {
    counts: Vec<u32>,
    freqs: Vec<f32>,
    n: u32,
}

impl<S> Recorder<S> for Frequency
where
    S: CountableState,
{
    fn initialize(state: &S) -> Self {
        Self {
            counts: vec![0; state.n_all()],
            freqs: vec![0.0; state.n_all()],
            n: 0,
        }
    }
    fn record(&mut self, state: &S) {
        for i in 0..state.n_all() {
            if state.get(i) {
                self.counts[i] += 1;
            }
        }
        self.n += 1;
    }

    fn finalize(&mut self) {
        self.freqs = self
            .counts
            .iter()
            .map(|&x| x as f32 / self.n as f32)
            .collect();
    }
}

pub trait Algorithm<M>
where
    M: Model,
{
    fn sample<P: Proposer<M::State>, R: Recorder<M::State>>(
        &mut self,
        state: &mut M::State,
        obs: &M::Observation,
        proposer: P,
    ) -> R;
}

pub struct MetropolisHasting<M> {
    model: M,
    iterations: usize,
    burn_in: usize,
}

impl<M> MetropolisHasting<M> {
    pub fn new(model: M, iterations: usize, burn_in: usize) -> Self {
        Self {
            model,
            iterations,
            burn_in,
        }
    }
}

impl<M> Algorithm<M> for MetropolisHasting<M>
where
    M: Model,
{
    fn sample<P, R>(&mut self, state: &mut M::State, obs: &M::Observation, proposer: P) -> R
    where
        P : Proposer<M::State>,
        R : Recorder<M::State>
    {
        let mut result = R::initialize(&state);

        let mut rng = rand::rng();

        for i in 0..(self.burn_in + self.iterations) {
            let m = proposer.propose(state, &mut rng);
            let log_q_ratio = proposer.log_proposal_ratio(state, &m);

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

            let log_l_ratio = match self.model.log_likelihood_ratio(state, obs, &m) {
                Some(log_p_ratio) => log_p_ratio,
                None => {
                    let log_l1 = self.model.log_likelihood(state, obs);
                    state.apply(&m);
                    let log_l2 = self.model.log_likelihood(state, obs);
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
