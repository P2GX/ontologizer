use crate::bayesian::state::{MgsaMove, MgsaState, State, ToggleSwap};
use indexmap::IndexSet;

/// A trait that connects a *State* configuration with *Observations* by assigning probabilities.
///
/// It defines the logic for the Prior P(S) and Likelihood P(O|S).
pub trait Model {
    type State: State;
    type Cache;

    /// Calculates the log prior probability P(S) of a state configuration.
    fn log_prior(&self, state: &Self::State) -> f64;

    /// Fast implementation for the log prior ratio P(S')/P(S).
    ///
    /// Returns `None` if an optimized calculation is not available for the specific move.
    fn log_prior_ratio(
        &self,
        _state: &Self::State,
        _m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }

    /// Calculates the log likelihood P(O | S) of the observations given the state.
    fn log_likelihood(&self, state: &Self::State, cache: &Self::Cache) -> f64;

    /// Fast implementation for the log likelihood ratio P(O|S')/P(O|S).
    ///
    /// Returns `None` if an optimized calculation is not available for the specific move.
    fn log_likelihood_ratio(
        &self,
        _state: &Self::State,
        _cache: &Self::Cache,
        _m: &<Self::State as State>::Move,
    ) -> Option<f64> {
        None
    }

    /// Creates a new cache based on the current state.
    fn create_cache(&self, state: &Self::State) -> Self::Cache;

    /// Updates the cache incrementally after a state change (move).
    fn update_cache(
        &self,
        cache: &mut Self::Cache,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    );

    /// Reverts the cache incrementally after reversing a state change.

    fn revert_cache(
        &self,
        cache: &mut Self::Cache,
        state: &Self::State,
        m: &<Self::State as State>::Move,
    );
}

#[derive(Debug, Clone)]
pub struct OrCache {
    /// Number of active terms targeting each gene (the latent signal).
    pub latent: Vec<usize>,

    /// Confusion Matrix Counts
    pub n_tp: usize, // Model = 1 | Observed = 1
    pub n_fp: usize, // Model = 0 | Observed = 1 (alpha)
    pub n_tn: usize, // Model = 0 | Observed = 0
    pub n_fn: usize, // Model = 1 | Observed = 0 (beta)
}

/// The "Noisy-OR" Bayesian Network model.
pub struct OrModel {
    terms_to_genes: Vec<IndexSet<usize>>, // Adjacency list (Term -> Genes)
    obs_genes: Vec<bool>,
    p_prior: (f64, f64),
    alpha_prior: (f64, f64),
    beta_prior: (f64, f64),
}

impl OrModel {
    pub fn new(
        terms_to_genes: &[IndexSet<usize>],
        obs_genes: &[bool],
        p: f64,
        alpha: f64,
        beta: f64,
    ) -> OrModel {
        // Define a "tightness" factor.
        // 1.0 = maximal variance for unimodal distribution (flattest curve possible).
        // 0.1 = sharp peak (10% of the max allowed variance).
        const VARIANCE_SCALE: f64 = 0.5;

        // Calculate safe priors that are guaranteed to be unimodal
        let p_prior = Self::set_unimodal_beta_prior(p, VARIANCE_SCALE);
        let alpha_prior = Self::set_unimodal_beta_prior(alpha, VARIANCE_SCALE);
        let beta_prior = Self::set_unimodal_beta_prior(beta, VARIANCE_SCALE);

        Self {
            terms_to_genes: terms_to_genes.to_vec(),
            obs_genes: obs_genes.to_vec(),
            p_prior,
            alpha_prior,
            beta_prior,
        }
    }

    /// Sets Beta parameters given a mean and a variance scaling factor.
    /// scale = 1.0 sets variance to the maximum possible value that preserves unimodality.
    /// scale < 1.0 sharper distributions.
    fn set_unimodal_beta_prior(mean: f64, scale: f64) -> (f64, f64) {
        // 1. Calculate the strict upper bound for variance to ensure alpha (a) > 1 and beta (b) > 1
        let max_unimodal_var = if mean <= 0.5 {
            (mean.powi(2) * (1.0 - mean)) / (1.0 + mean)
        } else {
            (mean * (1.0 - mean).powi(2)) / (2.0 - mean)
        };

        // 2. Apply the scaling factor (e.g., 50% of the limit)
        let var = max_unimodal_var * scale;
        // 3. Calculate shape parameters alpha (a), beta (b)
        // nu = [mu(1-mu) / var] - 1
        let nu = (mean * (1.0 - mean) / var) - 1.0;
        let a = mean * nu;
        let b = (1.0 - mean) * nu;

        (a, b)
    }

    #[cfg(test)]
    pub fn get_p_prior_params(&self) -> (f64, f64) {
        self.p_prior
    }
    #[cfg(test)]
    pub fn get_alpha_prior_params(&self) -> (f64, f64) {
        self.alpha_prior
    }
    #[cfg(test)]
    pub fn get_beta_prior_params(&self) -> (f64, f64) {
        self.beta_prior
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
                        if self.obs_genes[g] {
                            // Observed = 1
                            // FP -> TP
                            cache.n_fp -= 1;
                            cache.n_tp += 1;
                        } else {
                            // Observed = 0
                            // TN -> FN
                            cache.n_tn -= 1;
                            cache.n_fn += 1;
                        }
                    }
                    cache.latent[g] += 1;
                } else {
                    // Turning OFF: Latent count decreases
                    // Critical transition: 1 -> 0 (Gene becomes Predicted False)
                    if cache.latent[g] == 1 {
                        // Model 1 -> 0
                        if self.obs_genes[g] {
                            // Observed = 1
                            // TP -> FP
                            cache.n_tp -= 1;
                            cache.n_fp += 1;
                        } else {
                            // Observed = 0
                            // FN -> TN
                            cache.n_fn -= 1;
                            cache.n_tn += 1;
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

    // Log probability P(T|p)P(p)P(alpha)P(beta)
    fn log_prior(&self, state: &MgsaState) -> f64 {
        let n_active = state.terms.n_active() as f64;
        let n_inactive = state.terms.n_inactive() as f64;

        let p = state.params.p();
        let alpha = state.params.alpha();
        let beta = state.params.beta();

        let (p_a, p_b) = self.p_prior;
        let (alpha_a, alpha_b) = self.alpha_prior;
        let (beta_a, beta_b) = self.beta_prior;

        let log_prior_t = n_active * p.ln() + n_inactive * (1. - p).ln();
        let log_prior_p = (p_a - 1.0) * p.ln() + (p_b - 1.0) * (1.0 - p).ln();
        let log_prior_alpha = (alpha_a - 1.0) * alpha.ln() + (alpha_b - 1.0) * (1.0 - alpha).ln();
        let log_prior_beta = (beta_a - 1.0) * beta.ln() + (beta_b - 1.0) * (1.0 - beta).ln();
        log_prior_t + log_prior_p + log_prior_alpha + log_prior_beta
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
            MgsaMove::Parameter(inc) => match inc.index {
                0 => {
                    let n_active = state.terms.n_active() as f64;
                    let n_inactive = state.terms.n_inactive() as f64;

                    let p = state.params.p();
                    let p_new = p + inc.delta;

                    let (p_a, p_b) = self.p_prior;

                    let log_prior_ratio = n_active * (p_new / p).ln()
                        + n_inactive * ((1. - p_new) / (1. - p)).ln()
                        + (p_a - 1.0) * (p_new / p).ln()
                        + (p_b - 1.0) * ((1. - p_new) / (1. - p)).ln();
                    Some(log_prior_ratio)
                }
                1 => {
                    let alpha = state.params.alpha();
                    let alpha_new = state.params.alpha() + inc.delta;
                    let (alpha_a, alpha_b) = self.alpha_prior;

                    let log_prior_ratio = (alpha_a - 1.0) * (alpha_new / alpha).ln()
                        + (alpha_b - 1.0) * ((1. - alpha_new) / (1. - alpha)).ln();
                    Some(log_prior_ratio)
                }
                2 => {
                    let beta = state.params.beta();
                    let beta_new = state.params.beta() + inc.delta;
                    let (beta_a, beta_b) = self.beta_prior;

                    let log_prior_ratio = (beta_a - 1.0) * (beta_new / beta).ln()
                        + (beta_b - 1.0) * ((1. - beta_new) / (1. - beta)).ln();
                    Some(log_prior_ratio)
                }
                _ => None,
            },
        }
    }

    // Log probability P(O | H, alpha, beta).
    fn log_likelihood(&self, state: &MgsaState, cache: &OrCache) -> f64 {
        (cache.n_fp as f64) * state.params.alpha().ln()
            + (cache.n_tn as f64) * (1. - state.params.alpha()).ln()
            + (cache.n_fn as f64) * state.params.beta().ln()
            + (cache.n_tp as f64) * (1. - state.params.beta()).ln()
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

                let alpha = state.params.alpha();
                let beta = state.params.beta();

                let alpha_new = alpha + inc.delta;
                let beta_new = beta + inc.delta;

                let n_tp = cache.n_tp as f64;
                let n_fp = cache.n_fp as f64;
                let n_tn = cache.n_tn as f64;
                let n_fn = cache.n_fn as f64;

                let mut log_ll_ratio = 0.0;
                if inc.index == 1 {
                    // Only Alpha terms
                    log_ll_ratio += n_fp * (alpha_new / alpha).ln();
                    log_ll_ratio += n_tn * ((1.0 - alpha_new) / (1.0 - alpha)).ln();
                } else {
                    // Only Beta terms
                    log_ll_ratio += n_fn * (beta_new / beta).ln();
                    log_ll_ratio += n_tp * ((1.0 - beta_new) / (1.0 - beta)).ln();
                }

                Some(log_ll_ratio)
            }
        }
    }

    fn create_cache(&self, state: &MgsaState) -> OrCache {
        let n_genes = self.obs_genes.len();
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

        for (i, &obs) in self.obs_genes.iter().enumerate() {
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
            n_tp,
            n_fp,
            n_tn,
            n_fn,
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Test that the prior shaping logic correctly produces Beta distributions
    /// where both alpha (a) > 1 and beta (b) > 1, ensuring unimodality (bell curve shape).
    #[test]
    fn test_unimodal_prior_construction() {
        // Case 1: Mean 0.1 (Small)
        let (a, b) = OrModel::set_unimodal_beta_prior(0.1, 1.0);
        assert!(a > 1.0, "Alpha parameter must be > 1 for unimodality");
        assert!(b > 1.0, "Beta parameter must be > 1 for unimodality");

        let calculated_mean = a / (a + b);
        assert!((calculated_mean - 0.1).abs() < 1e-6);

        // Case 2: Mean 0.9 (Large)
        let (a, b) = OrModel::set_unimodal_beta_prior(0.9, 1.0);
        assert!(a > 1.0);
        assert!(b > 1.0);

        let calculated_mean = a / (a + b);
        assert!((calculated_mean - 0.9).abs() < 1e-6);

        // Case 3: Mean 0.5 (Middle)
        let (a, b) = OrModel::set_unimodal_beta_prior(0.5, 0.5); // Tighter variance
        assert!(a > 1.0);
        assert!(b > 1.0);
        assert_eq!(a, b, "For mean 0.5, a and b should be equal");
    }
}
