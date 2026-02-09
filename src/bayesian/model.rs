use crate::bayesian::proposer::{MgsaMove, ToggleSwap};
use crate::bayesian::state::{MgsaState, State};
use indexmap::IndexSet;
use rand::{Rng, RngExt};

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model {
    type State: State;
    type Cache;

    // Initialize the State
    fn initialize_state<R: Rng>(&self, rng: &mut R) -> Self::State;

    // Log probability P(S) to find a State configuration
    fn log_prior(&self, state: &Self::State) -> f64;

    // Fast implementation for log probability ratio P(S2)/P(S1). May not be available for the specific move.
    fn log_prior_ratio(
        &self,
        _state: &Self::State,
        _m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }

    fn create_cache(&self, state: &Self::State) -> Self::Cache;

    fn update_cache(
        &self,
        cache: &mut Self::Cache,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    );

    fn revert_cache(
        &self,
        cache: &mut Self::Cache,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    );

    // Log probability P(O | S) to find an Observable configuration given a State configuration.
    fn log_likelihood(&self, state: &Self::State, cache: &Self::Cache) -> f64;

    // Fast implementation for log likelihood ratio P(O|S2)/P(O|S1). May not be available for the specific move.
    fn log_likelihood_ratio(
        &self,
        state: &Self::State,
        _cache: &Self::Cache,
        _m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }
}

pub struct OrCache {
    pub latent: Vec<usize>,

    // Confusion Matrix Counts
    pub n_true_pos: usize,  // Model = 1 | Observed = 1
    pub n_false_pos: usize, // Model = 0 | Observed = 1 (alpha)
    pub n_true_neg: usize,  // Model = 0 | Observed = 0
    pub n_false_neg: usize, // Model = 1 | Observed = 0 (beta)
}

pub struct OrModel {
    p_init: f64,     // probability of term being active
    alpha_init: f64, // probability of a gene being incorrectly observed active.
    beta_init: f64,  // probability of a gene being incorrectly observed inactive.

    terms_to_genes: Vec<IndexSet<usize>>, // Adjacency list (Term -> Genes)
    observations: Vec<bool>,
}

impl OrModel {
    pub fn new(
        terms_to_genes: Vec<IndexSet<usize>>,
        observations: Vec<bool>,
        p: f64,
        alpha: f64,
        beta: f64,
    ) -> OrModel {
        Self {
            p_init: p,
            alpha_init: alpha,
            beta_init: beta,
            terms_to_genes,
            observations,
        }
    }

    /// Helper: Updates latent counts and confusion matrix based on a term toggle.
    fn update_cache_for_toggle(&self, cache: &mut OrCache, term_idx: usize, is_enabled: bool) {
        if let Some(genes) = self.terms_to_genes.get(term_idx) {
            for &g in genes {
                if is_enabled {
                    // Turning ON: Latent count increases
                    // n_tp: Model = 1 | Observed = 1
                    // n_fp: Model = 0 | Observed = 1
                    // n_tn: Model = 0 | Observed = 0
                    // n_fn: Model = 1 | Observed = 0

                    // Critical transition: 0 -> 1 (Gene becomes Predicted True)
                    if cache.latent[g] == 0 {
                        // Model 0 -> 1
                        if self.observations[g] {
                            // Observed = 1
                            // FP -> TP
                            cache.n_false_pos -= 1;
                            cache.n_true_pos += 1;
                        } else {
                            // Observed = 0
                            // TN -> FN
                            cache.n_true_neg -= 1;
                            cache.n_false_neg += 1;
                        }
                    }
                    cache.latent[g] += 1;
                } else {
                    // Turning OFF: Latent count decreases
                    // Critical transition: 1 -> 0 (Gene becomes Predicted False)
                    if cache.latent[g] == 1 {
                        // Model 1 -> 0
                        if self.observations[g] {
                            // Observed = 1
                            // TP -> FP
                            cache.n_true_pos -= 1;
                            cache.n_false_pos += 1;
                        } else {
                            // Observed = 0
                            // FN -> TN
                            cache.n_false_neg -= 1;
                            cache.n_true_neg += 1;
                        }
                    }
                    if cache.latent[g] > 0 {
                        cache.latent[g] -= 1;
                    } else {
                        panic!(
                            "Latent underflow detected for gene {}! State/Cache desync.",
                            g
                        );
                    }
                }
            }
        }
    }
}

impl Model for OrModel {
    type State = MgsaState;
    type Cache = OrCache;

    fn initialize_state<R: Rng>(&self, rng: &mut R) -> MgsaState {
        let n_terms = self.terms_to_genes.len();
        // Sample random terms based on fixed p
        let terms = (0..n_terms).map(|_| rng.random_bool(self.p_init)).collect();
        MgsaState::new(terms, self.p_init, self.alpha_init, self.beta_init)
    }

    // Log probability P(T|p)P(p)P(alpha)P(beta) to find a Term configuration
    fn log_prior(&self, state: &MgsaState) -> f64 {
        let m0 = state.terms.n_active() as f64;
        let m1 = state.terms.n_inactive() as f64;
        let p = state.params.p();

        let log_prior = m0 * p.ln() + m1 * (1. - p).ln();
        log_prior
    }

    fn log_prior_ratio(&self, state: &MgsaState, m: &MgsaMove) -> Option<f64> {
        match m {
            MgsaMove::Term(ts) => match ts {
                // log [P(T_new | p) / P(T_old | p)]
                ToggleSwap::Toggle(i) => {
                    let p = state.params.p();
                    if state.terms.get(*i) {
                        // Active -> Inactive
                        Some(((1.0 - p) / p).ln())
                    } else {
                        // Inactive -> Active
                        Some((p / (1.0 - p)).ln())
                    }
                }
                ToggleSwap::Swap(_, _) => Some(0.0),
            },
            MgsaMove::Parameter(inc) => {
                if inc.index == 0 {
                    // log [ P(T|p_new) / P(T|p_old)]
                    let p_old = state.params.p();
                    let p_new = p_old + inc.delta;

                    let m0 = state.terms.n_active() as f64;
                    let m1 = state.terms.n_inactive() as f64;

                    let log_prior_ratio =
                        m0 * (p_new / p_old).ln() + m1 * ((1. - p_new) / (1 - p_old)).ln();
                    Some(log_prior_ratio)
                } else {
                    Some(0.0)
                }
            }
        }
    }

    // Log probability P(O | H, alpha, beta)p(H|T) to find an observed Gene configuration given a Terms configuration.
    fn log_likelihood(&self, state: &MgsaState, cache: &OrCache) -> f64 {
        (cache.n_true_pos as f64) * (1. - state.params.beta()).ln()
            + (cache.n_false_pos as f64) * state.params.alpha().ln()
            + (cache.n_true_neg as f64) * (1. - state.params.alpha()).ln()
            + (cache.n_false_neg as f64) * state.params.beta().ln()
    }

    fn log_likelihood_ratio(
        &self,
        state: &Self::State,
        cache: &Self::Cache,
        m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        match m {
            MgsaMove::Term(_) => None,
            MgsaMove::Parameter(inc) => {
                // p changes don't affect Likelihood P(O|...)
                if inc.index == 0 {
                    return Some(0.0);
                }

                let alpha_old = state.params.alpha();
                let beta_old = state.params.beta();

                let mut alpha_new = alpha_old;
                let mut beta_new = beta_old;

                let n10 = cache.n_false_pos as f64;
                let n00 = cache.n_true_neg as f64;
                let n01 = cache.n_false_neg as f64;
                let n11 = cache.n_true_pos as f64;

                let mut log_ll_ratio = 0.0;
                if inc.index == 1 {
                    // Only Alpha terms
                    log_ll_ratio += n10 * (alpha_new / alpha_old).ln();
                    log_ll_ratio += n00 * ((1.0 - alpha_new) / (1.0 - alpha_old)).ln();
                } else {
                    // Only Beta terms
                    log_ll_ratio += n01 * (beta_new / beta_old).ln();
                    log_ll_ratio += n11 * ((1.0 - beta_new) / (1.0 - beta_old)).ln();
                }

                Some(log_ll_ratio)
            }
        }
    }

    fn create_cache(&self, state: &MgsaState) -> OrCache {
        let n_genes = self.observations.len();
        let mut latent = vec![0; n_genes];

        // Calculate Latent (O(N) full scan)
        // We need to iterate all active terms in the state
        for k in 0..state.terms.n_active() {
            let term_idx = state.terms.get_active(k);
            for &gene_idx in &self.terms_to_genes[term_idx] {
                latent[gene_idx] += 1;
            }
        }

        // Calculate Confusion Matrix
        let mut n_tp = 0;
        let mut n_fp = 0;
        let mut n_tn = 0;
        let mut n_fn = 0;

        for (i, &obs) in self.observations.iter().enumerate() {
            let predicted = latent[i] > 0;
            match (predicted, obs) {
                (true, true) => n_tp += 1,
                (false, true) => n_fp += 1,
                (true, false) => n_fn += 1,
                (false, false) => n_tn += 1,
            }
        }

        OrCache {
            latent,
            n_true_pos: n_tp,
            n_false_pos: n_fp,
            n_true_neg: n_tn,
            n_false_neg: n_fn,
        }
    }

    fn update_cache(&self, cache: &mut OrCache, updated_state: &MgsaState, m: &MgsaMove) {
        match m {
            MgsaMove::Term(ts) => match *ts {
                ToggleSwap::Toggle(i) => {
                    self.update_cache_for_toggle(cache, i, updated_state.terms.get(i));
                }
                ToggleSwap::Swap(i, j) => {
                    self.update_cache_for_toggle(cache, i, updated_state.terms.get(i));
                    self.update_cache_for_toggle(cache, j, updated_state.terms.get(j));
                }
            },
            MgsaMove::Parameter(_) => {
                // Parameters (p, alpha, beta) do not affect latent counts or confusion matrix structure.
            }
        }
    }

    fn revert_cache(&self, cache: &mut OrCache, state: &MgsaState, m: &MgsaMove) {
        self.update_cache(cache, state, m);
    }
}
