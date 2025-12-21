use crate::bayesian::proposer::ToggleSwap;
use crate::bayesian::proposer::ToggleSwap::{Swap, Toggle};
use crate::bayesian::state::{BinaryParameterState, State};
use rand::Rng;

// A trait that connects *STATE* and *OBSERVATION* by assigning probabilities.
pub trait Model {
    type State: State;
    type Cache;

    // Initialize the State
    fn sample_prior<R: Rng>(&self, rng: &mut R, n: usize) -> Vec<<Self::State as State>::Value>;

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
    fn log_likelihood(&self, cache: &Self::Cache) -> f64;

    // Fast implementation for log likelihood ratio P(O|S2)/P(O|S1). May not be available for the specific move.
    fn log_likelihood_ratio(
        &self,
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

pub struct OrModel<S> {
    p: f64,     // probability of term being active
    alpha: f64, // probability of a gene being incorrectly observed active.
    beta: f64,  // probability of a gene being incorrectly observed inactive.

    terms_to_genes: Vec<Vec<usize>>, // Adjacency list (Term -> Genes)
    observations: Vec<bool>,
    _marker: std::marker::PhantomData<S>,
}

impl<S> OrModel<S> {
    pub fn new(
        terms_to_genes: Vec<Vec<usize>>,
        observations: Vec<bool>,
        p: f64,
        alpha: f64,
        beta: f64,
    ) -> OrModel<S> {
        Self {
            p,
            alpha,
            beta,
            terms_to_genes,
            observations,
            _marker: std::marker::PhantomData,
        }
    }

    fn update_toggle(&self, cache: &mut OrCache, term_idx: usize, enable: bool) {
        if let Some(genes) = self.terms_to_genes.get(term_idx) {
            for &g in genes {
                if enable {
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
                    cache.latent[g] -= 1;
                }
            }
        }
    }

    pub fn heuristic_start(&self) -> Vec<bool> {
        let mut prior = vec![false; self.terms_to_genes.len()];
        for (term_index, gene_indices) in self.terms_to_genes.iter().enumerate() {
            for &gene_index in gene_indices {
                if self.observations[gene_index] {
                    prior[term_index] = true;
                    break;
                }
            }
        }
        prior
    }
}

impl<S> Model for OrModel<S>
where
    S: BinaryParameterState<Move = ToggleSwap>,
{
    type State = S;
    type Cache = OrCache;

    fn sample_prior<R: Rng>(&self, rng: &mut R, n_terms: usize) -> Vec<bool> {
        (0..n_terms).map(|_| rng.random_bool(self.p)).collect()
    }

    // Log probability P(T) to find a Term configuration
    fn log_prior(&self, state: &S) -> f64 {
        let m0 = state.n_active() as f64;
        let m1 = state.n_inactive() as f64;

        m0 * self.p.ln() + m1 * (1. - self.p).ln()
    }

    fn log_prior_ratio(&self, state: &Self::State, m: &ToggleSwap) -> Option<f64> {
        let log_ratio = match *m {
            Toggle(i) => {
                if state.get(i) {
                    // from active to inactive
                    ((1.0 - self.p) / (self.p)).ln()
                } else {
                    // from inactive to active
                    (self.p / (1.0 - self.p)).ln()
                }
            }
            Swap(_, _) => 0.0,
        };
        Some(log_ratio)
    }

    fn create_cache(&self, state: &S) -> OrCache {
        let n_genes = self.observations.len();
        let mut latent = vec![0; n_genes];

        // Calculate Latent (O(N) full scan)
        // We need to iterate all active terms in the state
        for k in 0..state.n_active() {
            let term_idx = state.get_active(k);
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

    fn update_cache(&self, cache: &mut Self::Cache, state: &Self::State, m: &ToggleSwap) {
        match *m {
            Toggle(i) => {
                self.update_toggle(cache, i, state.get(i));
            }
            Swap(i, j) => {
                self.update_toggle(cache, i, state.get(i));
                self.update_toggle(cache, j, state.get(j));
            }
        }
    }

    fn revert_cache(&self, cache: &mut Self::Cache, state: &Self::State, m: &<Self::State as State>::Move) {
        self.update_cache(cache, state, m);
    }
    
    
    // Log probability P(O | T) to find an observed Gene configuration given a Terms configuration.
    fn log_likelihood(&self, cache: &OrCache) -> f64 {
        (cache.n_true_pos as f64) * (1. - self.beta).ln()
            + (cache.n_false_pos as f64) * self.alpha.ln()
            + (cache.n_true_neg as f64) * (1. - self.alpha).ln()
            + (cache.n_false_neg as f64) * self.beta.ln()
    }
}
