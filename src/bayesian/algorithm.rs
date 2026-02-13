use crate::bayesian::model::Model;
use crate::bayesian::proposer::Proposer;
use crate::bayesian::recorder::Recorder;
use crate::bayesian::state::State;
use rand::RngExt;

/// A generic interface for Bayesian sampling algorithms.
pub trait Algorithm<M>
where
    M: Model,
{
    fn sample<R: Recorder<M::State>>(&mut self, state: &mut M::State) -> R::Target;
}

pub struct MetropolisHastings<M: Model, P> {
    model: M,
    proposer: P,
    n_iterations: usize,
    burn_in: usize,
}

/// The Metropolis-Hastings algorithm for sampling from probability distributions.
///
/// It constructs a Markov chain that asymptotically converges to the desired distribution.
impl<M, P> MetropolisHastings<M, P>
where
    M: Model,
    P: Proposer<M::State>,
{
    pub fn new(model: M, proposer: P, n_iterations: usize, burn_in: usize) -> Self {
        Self {
            model,
            proposer,
            n_iterations,
            burn_in,
        }
    }

    /// Calculates the log proposal ratio: log( q(x'|x) / q(x|x') )
    fn get_log_proposal_ratio(
        &mut self,
        state: &mut M::State,
        m: &<M::State as State>::Move,
    ) -> f64 {
        match self.proposer.log_proposal_ratio(state, &m) {
            Some(log_q_ratio) => log_q_ratio,
            None => {
                let log_q_reverse = self.proposer.log_proposal(state);
                state.apply(&m);
                let log_q_forward = self.proposer.log_proposal(state);
                state.revert(&m);
                log_q_forward - log_q_reverse
            }
        }
    }

    /// Calculates the log prior ratio: log( P(x') / P(x) )
    fn get_log_prior_ratio(&mut self, state: &mut M::State, m: &<M::State as State>::Move) -> f64 {
        match self.model.log_prior_ratio(state, &m) {
            Some(log_p_ratio) => log_p_ratio,
            None => {
                let log_p_current = self.model.log_prior(state);
                state.apply(&m);
                let log_p_proposed = self.model.log_prior(state);
                state.revert(&m);
                log_p_proposed - log_p_current
            }
        }
    }

    /// Calculates the log likelihood ratio: log( L(x') / L(x) )
    fn get_log_likelihood_ratio(
        &mut self,
        state: &mut M::State,
        cache: &mut M::Cache,
        m: &<M::State as State>::Move,
    ) -> f64 {
        match self.model.log_likelihood_ratio(state, cache, m) {
            Some(log_l_ratio) => log_l_ratio,
            None => {
                let log_l_current = self.model.log_likelihood(state, cache);
                state.apply(m);
                self.model.update_cache(cache, state, m);
                let log_l_proposed = self.model.log_likelihood(state, cache);
                state.revert(m);
                self.model.revert_cache(cache, state, m);
                log_l_proposed - log_l_current
            }
        }
    }
}

impl<M, P> Algorithm<M> for MetropolisHastings<M, P>
where
    M: Model,
    P: Proposer<M::State>,
{
    fn sample<R>(&mut self, state: &mut M::State) -> R::Target
    where
        R: Recorder<M::State>,
    {
        let mut rng = rand::rng();

        // ------ TRANSIENT PHASE (BURN-IN) ------
        // Run the chain to reach the stationary distribution but do not record.
        let mut cache = self.model.create_cache(state);

        for _ in 0..self.burn_in {
            let m = self.proposer.propose(state, &mut rng);

            let log_q_ratio = self.get_log_proposal_ratio(state, &m);
            let log_p_ratio = self.get_log_prior_ratio(state, &m);
            let log_l_ratio = self.get_log_likelihood_ratio(state, &mut cache, &m);

            // Acceptance probability alpha = min(1, likelihood * prior / proposal)
            // log_alpha = log_l + log_p - log_q
            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let random: f64 = rng.random_range(0.0..1.0);
            let log_random: f64 = random.ln();
            if log_random < log_accept {
                state.apply(&m);
                self.model.update_cache(&mut cache, state, &m);
            }
        }

        // ------ STATIONARY PHASE (SAMPLING) ------
        // Run the chain and record.
        let mut result = R::initialize(&state);

        for i in 0..self.n_iterations {
            let m = self.proposer.propose(state, &mut rng);

            let log_q_ratio = self.get_log_proposal_ratio(state, &m);
            let log_p_ratio = self.get_log_prior_ratio(state, &m);
            let log_l_ratio = self.get_log_likelihood_ratio(state, &mut cache, &m);

            let log_accept = log_l_ratio + log_p_ratio - log_q_ratio;

            let random: f64 = rng.random_range(0.0..1.0);
            let log_random: f64 = random.ln();
            if log_random < log_accept {
                state.apply(&m);
                self.model.update_cache(&mut cache, state, &m);
                result.record(&m, i);
            }
        }
        result.finalize(self.n_iterations)
    }
}

// ==========================================
// TESTS
// ==========================================
#[cfg(test)]
mod test_parameter_inference {
    use super::*;
    use crate::bayesian::measure::Mean;
    use crate::bayesian::model::{Model, OrModel};
    use crate::bayesian::proposer::ParameterGaussProposer;
    use crate::bayesian::recorder::ParameterRecorder;
    use crate::bayesian::state::MgsaState;
    use crate::core::result::Measure;
    use indexmap::{IndexSet, indexset};

    fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn test_parameter_inference(
        terms: Vec<bool>,
        obs_genes: &Vec<bool>,
        terms_to_obs_genes: &Vec<IndexSet<usize>>,
        p: f64,
        alpha: f64,
        beta: f64,
    ) {
        let model = OrModel::new(
            terms_to_obs_genes.clone(),
            obs_genes.clone(),
            p,
            alpha,
            beta,
        );
        let mut state = MgsaState::new(terms, p, alpha, beta);

        // Calculate theoretical posteriors based on counts
        let cache = model.create_cache(&state);
        let p_prior = model.get_p_prior_params();
        let alpha_prior = model.get_alpha_prior_params();
        let beta_prior = model.get_beta_prior_params();

        let n_active = state.n_terms_active();
        let n_inactive = state.n_terms_inactive();

        let p_post = (p_prior.0 + n_active as f64, p_prior.1 + n_inactive as f64);
        let alpha_post = (
            alpha_prior.0 + cache.n_fp as f64,
            alpha_prior.1 + cache.n_tn as f64,
        );
        let beta_post = (
            beta_prior.0 + cache.n_fn as f64,
            beta_prior.1 + cache.n_tp as f64,
        );

        let post = vec![
            p_post.0 / (p_post.0 + p_post.1),
            alpha_post.0 / (alpha_post.0 + alpha_post.1),
            beta_post.0 / (beta_post.0 + beta_post.1),
        ];

        let proposer = ParameterGaussProposer::new(1.0);
        let mut algorithm = MetropolisHastings::new(model, proposer, 500_000, 50_000);

        let measure: Vec<Mean> = algorithm.sample::<ParameterRecorder>(&mut state);

        println!(
            "Measured: {:?}",
            measure.iter().map(|m| m.score()).collect::<Vec<_>>()
        );
        println!("Theoretical: {:?}", post);

        for (measure, theo) in measure.iter().zip(post) {
            assert!(approx_equal(measure.score(), theo, 0.05));
        }
    }

    #[test]
    fn test_parameter_small() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let terms = vec![false, false, false, true];
        let terms_to_genes: Vec<IndexSet<usize>> =
            vec![indexset! {0}, indexset! {1}, indexset! {2}, indexset! {3}];
        // Observed gene state
        let obs_genes = vec![false, false, false, true];

        test_parameter_inference(terms, &obs_genes, &terms_to_genes, p, alpha, beta);
    }

    #[test]
    fn test_parameter_large() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let terms = vec![true, true, false, false, false, false, false, false, false];

        // Term 8 activates genes 0, 1, 2
        // Term 9 activates genes 3, 4
        // (Total Hidden Active = 5)
        let terms_to_genes: Vec<IndexSet<usize>> = vec![
            indexset! {0, 1, 2}, // Active Term 1
            indexset! {3, 4},    // Active Term 2
            indexset! {5, 6},    // Targets false positive region
            indexset! {7},
            indexset! {8, 9},
            indexset! {10, 11}, // Inactive term noise
            indexset! {12, 13},
            indexset! {14, 15},
            indexset! {18, 19},
            indexset! {16, 17},
        ];

        // Observed State: Constructed to match counts
        // TP=3, FN=2, FP=3, TN=12
        let obs_genes = vec![
            true, true, true, // 0-2: True Positives (3)
            false, false, // 3-4: False Negatives (2)
            true, true, true, // 5-7: False Positives (3)
            false, false, // 8-9: True Negatives (2)
            false, false, false, false, false, // 10-14: True Negatives (5)
            false, false, false, false, false, // 15-19: True Negatives (5)
        ];

        test_parameter_inference(terms, &obs_genes, &terms_to_genes, p, alpha, beta);
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use crate::bayesian::measure::Probability;
    use crate::bayesian::model::OrModel;
    use crate::bayesian::proposer::{TermToggleProposer, TermToggleSwapProposer};
    use crate::bayesian::recorder::TermRecorder;
    use crate::bayesian::state::MgsaState;
    use crate::core::result::Measure;
    use indexmap::{IndexSet, indexset};

    fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn test_term_inference(
        terms: Vec<bool>,
        obs_genes: &Vec<bool>,
        terms_to_genes: &Vec<IndexSet<usize>>,
        p: f64,
        alpha: f64,
        beta: f64,
        posterior: Vec<f64>,
    ) {
        let model = OrModel::new(terms_to_genes.clone(), obs_genes.clone(), p, alpha, beta);
        let mut state = MgsaState::new(terms, p, alpha, beta);
        let proposer = TermToggleProposer::new();
        let mut algorithm = MetropolisHastings::new(model, proposer, 10_000, 0);

        let measure: Vec<Probability> = algorithm.sample::<TermRecorder>(&mut state);

        println!(
            "Measured: {:?}",
            measure.iter().map(|m| m.score()).collect::<Vec<_>>()
        );
        println!("Theoretical: {:?}", posterior);

        for (measure, theo) in measure.iter().zip(posterior) {
            assert!(approx_equal(measure.score(), theo, 0.05));
        }
    }

    #[test]
    fn test_one_term_one_gene() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let terms = vec![false];
        let terms_to_genes: Vec<IndexSet<usize>> = vec![indexset! {0}];

        // Active Observation O=(1)
        let obs_genes = vec![true];
        let post = vec![p * (1. - beta) / (alpha * (1. - p) + (1. - beta) * p)];

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            post,
        );

        // Inactive Observation O=(0)
        let obs_genes = vec![false];
        let post = vec![p * beta / ((1. - alpha) * (1. - p) + beta * p)];

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            post,
        );
    }

    #[test]
    fn test_two_terms_one_gene() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let terms = vec![false, false];
        let terms_to_genes = vec![indexset! {0}, indexset! {0}];

        // Case 1: Active Observation O=(1)
        let obs_genes = vec![true];
        // Calculation for posterior...
        let numerator = p * (1. - beta);
        let denominator = alpha * (1. - p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.));
        let posterior_val = numerator / denominator;

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            vec![posterior_val, posterior_val],
        );

        // Case 2: Inactive Observation O=(0)
        let observations = vec![false];
        let numerator = p * beta;
        let denominator = (1. - alpha) * (1. - p).powf(2.) + beta * (1. - (1. - p).powf(2.));
        let posterior_val = numerator / denominator;

        test_term_inference(
            terms,
            &observations,
            &terms_to_genes,
            p,
            alpha,
            beta,
            vec![posterior_val, posterior_val],
        );
    }

    #[test]
    fn test_full_two_term_two_gene() {
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let terms = vec![false, false];
        let terms_to_genes: Vec<IndexSet<usize>> = vec![indexset! {0, 1}, indexset! {0, 1}];

        // Case 1: O=(1, 1)
        let obs_genes = vec![true, true];
        let num = p * (1. - beta).powf(2.);
        let den =
            alpha.powf(2.) * (1. - p).powf(2.) + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.));
        let post = num / den;

        let posterior = vec![post, post];

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            posterior,
        );

        // Case 2: O=(0, 0)
        let obs_genes = vec![false, false];
        let num = p * beta.powf(2.);
        let den =
            (1. - alpha).powf(2.) * (1. - p).powf(2.) + beta.powf(2.) * (1. - (1. - p).powf(2.));
        let post = num / den;

        let posterior = vec![post, post];

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
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

        let terms = vec![false, false];
        let terms_to_genes: Vec<IndexSet<usize>> = vec![
            indexset! {0},    // Term 0
            indexset! {0, 1}, // Term 1
        ];

        // Case 1: O=(1, 0)
        let obs_genes = vec![true, false];
        // Likelihoods for states (00, 10, 01, 11)
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

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );

        // Case 2: O=(0, 1)
        let obs_genes = vec![false, true];
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

        test_term_inference(
            terms.clone(),
            &obs_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );
    }

    #[test]
    fn test_term_broadness() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let terms = vec![false, false, false];
        let terms_to_genes: Vec<IndexSet<usize>> = vec![
            indexset! {0, 1, 2},
            indexset! {0, 1, 2, 3, 4},
            indexset! {0, 1, 2, 3, 4, 5, 6},
        ];

        let obs_genes = vec![true, true, true, false, false, false, false];

        let model = OrModel::new(terms_to_genes.clone(), obs_genes.clone(), p, alpha, beta);
        let mut state = MgsaState::new(terms, p, alpha, beta);
        let proposer = TermToggleProposer::new();
        let mut algorithm = MetropolisHastings::new(model, proposer, 10_000, 0);

        let measure: Vec<Probability> = algorithm.sample::<TermRecorder>(&mut state);

        // Ensure specific terms are more likely than others based on specificity
        assert!(measure[0].score() > measure[1].score());
        assert!(measure[1].score() > measure[2].score());
    }

    #[test]
    fn test_term_background_noise() {
        fn log_sum_exp(vals: &[f64]) -> f64 {
            let max_val = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let sum = vals.iter().map(|v| (v - max_val).exp()).sum::<f64>();
            max_val + sum.ln()
        }

        let n = 500;
        let p = 0.01;
        let alpha = 0.05;
        let beta = 0.10;

        let terms = vec![false; n];
        let mut terms_to_genes: Vec<IndexSet<usize>> = (0..n).map(|_| (0..n).collect()).collect();
        let mut obs_genes = vec![false; n];

        // Setup specific case: gene 0 is true, term 0 points only to gene 0
        obs_genes[0] = true;
        terms_to_genes[0] = indexset! {0};

        let model = OrModel::new(terms_to_genes.clone(), obs_genes.clone(), p, alpha, beta);
        let mut state = MgsaState::new(terms, p, alpha, beta);
        let proposer = TermToggleSwapProposer::new();
        let mut algorithm = MetropolisHastings::new(model, proposer, 500_000, 100_000);

        let measure: Vec<Probability> = algorithm.sample::<TermRecorder>(&mut state);

        // --- Manual Posterior Calculation for Verification ---
        let n_background = n - 1; // 49
        let ln = |x: f64| x.ln();

        // Priors (Log space)
        // P(t1=1) and P(t1=0)
        let lp_t1_1 = ln(p);
        let lp_t1_0 = ln(1.0 - p);

        // P(S=0) = (1-p)^499
        let lp_s_0 = (n_background as f64) * ln(1.0 - p);
        // P(S=1) = 1 - (1-p)^499. computed in linear space for precision then logged
        let p_s_1 = 1.0 - (1.0 - p).powi(n_background as i32);
        let lp_s_1 = ln(p_s_1);

        let ln_alpha = ln(alpha);
        let ln_not_alpha = ln(1.0 - alpha);
        let ln_beta = ln(beta);
        let ln_not_beta = ln(1.0 - beta); // 1-beta is sensitivity (0.9)

        let ll_00 = ln_alpha + (n_background as f64) * ln_not_alpha;
        let ll_10 = ln_not_beta + (n_background as f64) * ln_not_alpha;

        // Likelihoods
        let ll_00 = ln_alpha + (n_background as f64) * ln_not_alpha;
        let ll_10 = ln_not_beta + (n_background as f64) * ln_not_alpha;
        let ll_01 = ln_not_beta + (n_background as f64) * ln_beta;
        let ll_11 = ll_01;

        // Joint
        let lj_00 = ll_00 + lp_t1_0 + lp_s_0;
        let lj_10 = ll_10 + lp_t1_1 + lp_s_0;
        let lj_01 = ll_01 + lp_t1_0 + lp_s_1;
        let lj_11 = ll_11 + lp_t1_1 + lp_s_1;

        let log_evidence = log_sum_exp(&[lj_00, lj_10, lj_01, lj_11]);

        let prob_t1 = (log_sum_exp(&[lj_10, lj_11]) - log_evidence).exp();

        // P(S=1 | o)
        let prob_s = (log_sum_exp(&[lj_01, lj_11]) - log_evidence).exp();
        // P(tk=1 | o) approximation
        let prob_tk = prob_s * (p / p_s_1);

        println!(
            "Term 0 Posterior: Theo={:.4}, Sim={:.4}",
            prob_t1,
            measure[0].score()
        );
        println!(
            "Background Term Posterior: Theo={:.4}, Sim={:.4}",
            prob_tk,
            measure[1].score()
        );

        assert!(approx_equal(prob_t1, measure[0].score(), 0.05));
        assert!(approx_equal(prob_tk, measure[1].score(), 0.05));
    }
}
