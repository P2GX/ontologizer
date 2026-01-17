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
    fn sample<R: Recorder<M::State>>(&mut self, state: &mut M::State) -> R::Target;
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
        cache: &mut M::Cache,
        m: &<M::State as State>::Move,
    ) -> f64 {
        match self.model.log_likelihood_ratio(cache, m) {
            Some(log_l_ratio) => log_l_ratio,
            None => {
                let log_l1 = self.model.log_likelihood(cache);
                state.apply(m);
                self.model.update_cache(cache, state, m);
                let log_l2 = self.model.log_likelihood(cache);
                state.revert(m);
                self.model.revert_cache(cache, state, m);
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
    fn sample<R>(&mut self, state: &mut M::State) -> R::Target
    where
        R: Recorder<M::State>,
    {
        let mut rng = rand::rng();

        // ------ TRANSIENT ------
        let mut cache = self.model.create_cache(state);
        for _ in 0..self.burn_in {
            let m = self.proposer.propose(state, &mut rng);
            let log_q_ratio = self.get_log_proposal_ratio(state, &m);
            let log_p_ratio = self.get_log_prior_ratio(state, &m);
            let log_l_ratio = self.get_log_likelihood_ratio(state, &mut cache, &m);
            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let x: f64 = rng.random_range(0.0..1.0);
            let accept = log_accept >= 0.0 || x.ln() < log_accept;
            if accept {
                state.apply(&m);
                self.model.update_cache(&mut cache, state, &m);
            }
        }

        // ------ (HOPEFULLY) STATIONARY ------
        let mut result = R::initialize(&state);
        for i in 0..self.iterations {
            let m = self.proposer.propose(state, &mut rng);
            let log_q_ratio = self.get_log_proposal_ratio(state, &m);
            let log_p_ratio = self.get_log_prior_ratio(state, &m);
            let log_l_ratio = self.get_log_likelihood_ratio(state, &mut cache, &m);
            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let x: f64 = rng.random_range(0.0..1.0);
            let accept = log_accept >= 0.0 || x.ln() < log_accept;
            if accept {
                state.apply(&m);
                self.model.update_cache(&mut cache, state, &m);
                result.record(&m, i);
            }
        }
        result.finalize(self.iterations)
    }
}

#[cfg(test)]
mod test {
    use crate::bayesian::measure::Probability;
    use super::*;
    use crate::bayesian::model::OrModel;
    use crate::bayesian::proposer::{UniformProposer, UniformToggleProposer};
    use crate::bayesian::recorder::ProbabilityRecorder;
    use crate::bayesian::state::MgsaState;
    use crate::core::result::Measure;

    fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn test_network(
        raw_state: Vec<bool>,
        observations: &Vec<bool>,
        state_to_observations: &Vec<Vec<usize>>,
        p: f64,
        alpha: f64,
        beta: f64,
        posterior: Vec<f64>,
    ) {
        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(raw_state);
        let proposer = UniformToggleProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Probability = algorithm.sample::<ProbabilityRecorder>(&mut state);
        println! {"{:?}", measure.scores().collect::<Vec<_>>()}
        for (sim, theo) in measure.scores().zip(posterior) {
            assert!(approx_equal(sim, theo, 0.05));
        }
    }

    #[test]
    fn test_background_noise() {
        fn log_sum_exp(vals: &[f64]) -> f64 {
            let max_val = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let sum = vals.iter().map(|v| (v - max_val).exp()).sum::<f64>();
            max_val + sum.ln()
        }

        let n = 500;

        let p = 0.01;
        let alpha = 0.05;
        let beta = 0.10;

        let state = vec![false; n];
        let mut state_to_observations: Vec<Vec<usize>> = (0..n).map(|_| (0..n).collect()).collect();

        // Observations
        let mut observations = vec![false; n];

        observations[0] = true;

        state_to_observations[0] = vec![0];

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 500_000, 100_000);

        let measure: Probability = algorithm.sample::<ProbabilityRecorder>(&mut state);

        let n_background = n - 1; // 49

        // Log-probabilities helper
        let ln = |x: f64| x.ln();

        // Priors (Log space)
        // P(t1=1) and P(t1=0)
        let lp_t1_1 = ln(p);
        let lp_t1_0 = ln(1.0 - p);

        // P(S=0) = (1-p)^49
        let lp_s_0 = (n_background as f64) * ln(1.0 - p);
        // P(S=1) = 1 - (1-p)^49. computed in linear space for precision then logged
        let p_s_1 = 1.0 - (1.0 - p).powi(n_background as i32);
        let lp_s_1 = ln(p_s_1);

        // Likelihoods (Log space)
        // Common factors
        let ln_alpha = ln(alpha);
        let ln_not_alpha = ln(1.0 - alpha);
        let ln_beta = ln(beta);
        let ln_not_beta = ln(1.0 - beta); // 1-beta is sensitivity (0.9)

        // L_00: t1=0, S=0. h1=0, h_rest=0.
        // o1=1 (FP), o_rest=0 (TN)
        let ll_00 = ln_alpha + (n_background as f64) * ln_not_alpha;

        // L_10: t1=1, S=0. h1=1, h_rest=0.
        // o1=1 (TP), o_rest=0 (TN)
        let ll_10 = ln_not_beta + (n_background as f64) * ln_not_alpha;

        // L_01: t1=0, S=1. h1=1, h_rest=1.
        // o1=1 (TP), o_rest=0 (FN)
        let ll_01 = ln_not_beta + (n_background as f64) * ln_beta;

        // L_11: t1=1, S=1. h1=1, h_rest=1.
        // Same likelihood as L_01
        let ll_11 = ll_01;

        // Joint Log Probabilities: Likelihood + Prior
        let lj_00 = ll_00 + lp_t1_0 + lp_s_0;
        let lj_10 = ll_10 + lp_t1_1 + lp_s_0;
        let lj_01 = ll_01 + lp_t1_0 + lp_s_1;
        let lj_11 = ll_11 + lp_t1_1 + lp_s_1;

        // Log-Sum-Exp to get Evidence E
        let log_evidence = log_sum_exp(&[lj_00, lj_10, lj_01, lj_11]);

        // Posteriors
        // P(t1=1 | o) = (J_10 + J_11) / E
        let log_prob_t1 = log_sum_exp(&[lj_10, lj_11]) - log_evidence;
        let prob_t1 = log_prob_t1.exp();

        // P(S=1 | o)
        let log_prob_s = log_sum_exp(&[lj_01, lj_11]) - log_evidence;
        let prob_s = log_prob_s.exp();

        // P(tk=1 | o) = P(S=1 | o) * p / P(S=1)
        let prob_tk = prob_s * (p / p_s_1);

        println!("{:?} {:?}", prob_t1, measure.get_score(0));
        println!("{:?} {:?}", prob_tk, measure.get_score(1));
        assert!(approx_equal(prob_t1, measure.get_score(0), 0.05));
        assert!(approx_equal(prob_tk, measure.get_score(1), 0.05));
    }

    #[test]
    fn test_term_broadness() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let raw_state = vec![false, false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0, 1, 2],
            vec![0, 1, 2, 3, 4],
            vec![0, 1, 2, 3, 4, 5, 6],
        ];

        // Observations
        let observations = vec![true, true, true, false, false, false, false];

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(raw_state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Probability = algorithm.sample::<ProbabilityRecorder>(&mut state);
        println! {"{:?}",measure.scores().collect::<Vec<_>>()}
        assert!(measure.get_score(0) > measure.get_score(1));
        assert!(measure.get_score(1) > measure.get_score(2));
    }

    #[test]
    fn test_one_term_one_gene() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let raw_state = vec![false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0]];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior = vec![p * (1. - beta) / (alpha * (1. - p) + (1. - beta) * p)];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false];
        let posterior = vec![p * beta / ((1. - alpha) * (1. - p) + beta * p)];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_two_term_one_gene() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0], vec![0]];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior_t1 =
            p * (1. - beta) / (alpha * (1. - p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));
        let posterior_t2 =
            p * (1. - beta) / (alpha * (1. - p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false];
        let posterior_t1 =
            p * beta / ((1. - alpha) * (1. - p).powf(2.) + beta * (1. - (1. - p).powf(2.)));
        let posterior_t2 =
            p * beta / ((1. - alpha) * (1. - p).powf(2.) + beta * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_full_two_term_two_gene() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0, 1], vec![0, 1]];

        // Active Observation O=(1)
        let observations = vec![true, true];
        let posterior_t1 = p * (1. - beta).powf(2.)
            / (alpha.powf(2.) * (1. - p).powf(2.)
            + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * (1. - beta).powf(2.)
            / (alpha.powf(2.) * (1. - p).powf(2.)
            + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false, false];
        let posterior_t1 = p * beta.powf(2.)
            / ((1. - alpha).powf(2.) * (1. - p).powf(2.)
            + beta.powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * beta.powf(2.)
            / ((1. - alpha).powf(2.) * (1. - p).powf(2.)
            + beta.powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_asymmetric_two_term_two_gene() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0],    // Term 0
            vec![0, 1], // Term 1
        ];

        let observations = vec![true, false];

        let lh_00 = alpha * (1.0 - alpha);
        let lh_10 = (1.0 - beta) * (1.0 - alpha);
        let lh_01 = (1.0 - beta) * beta;
        let lh_11 = (1.0 - beta) * beta;

        let joint_00 = lh_00 * (1.0 - p).powf(2.0);
        let joint_10 = lh_10 * p * (1.0 - p);
        let joint_01 = lh_01 * p * (1.0 - p);
        let joint_11 = lh_11 * p.powf(2.0);

        let z = joint_00 + joint_10 + joint_01 + joint_11;

        // Marginalize for Posteriors
        let post_t0 = (joint_10 + joint_11) / z;
        let post_t1 = (joint_01 + joint_11) / z;

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );

        let observations = vec![false, true];

        let lh_00 = (1.0 - alpha) * alpha;
        let lh_10 = beta * alpha;
        let lh_01 = beta * (1.0 - beta);
        let lh_11 = beta * (1.0 - beta);

        let joint_00 = lh_00 * (1.0 - p).powf(2.0);
        let joint_10 = lh_10 * p * (1.0 - p);
        let joint_01 = lh_01 * p * (1.0 - p);
        let joint_11 = lh_11 * p.powf(2.0);

        let z = joint_00 + joint_10 + joint_01 + joint_11;

        let post_t0 = (joint_10 + joint_11) / z;
        let post_t1 = (joint_01 + joint_11) / z;

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );
    }
}