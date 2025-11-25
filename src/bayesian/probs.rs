use crate::core::GeneSymbol;
use ontolius::TermId;
use std::collections::HashMap;

pub trait Probability {
    fn probability() -> f64;

    fn log_probability() -> f64;
}

pub fn neighborhood_size(terms: &HashMap<TermId, bool>) -> usize {
    let m = terms.len();
    let m0 = terms.iter().filter(|&(k, &v)| v == false).count();
    let m1 = terms.iter().filter(|&(k, &v)| v == true).count();
    m + (m0 * m1)
}

/// Compute n00, n01, n10, n11 counts for the likelihood.
/// Requires the global gene universe (e.g. all annotated genes).
pub fn counts_nxy(obs: &HashMap<GeneSymbol, bool>, hid: &HashMap<GeneSymbol, bool>) -> [u32; 4] {
    let (mut n00, mut n01, mut n10, mut n11) = (0, 0, 0, 0);
    for g in obs.keys() {
        let o = obs[g];
        let h = hid[g];
        match (o, h) {
            (false, false) => n00 += 1,
            (false, true) => n01 += 1,
            (true, false) => n10 += 1,
            (true, true) => n11 += 1,
        }
    }
    [n00, n01, n10, n11]
}

pub fn log_prior(q: f64, terms: &HashMap<TermId, bool>) -> f64 {
    let m0 = terms.iter().filter(|&(_, &v)| !v).count() as f64;
    let m1 = terms.iter().filter(|&(_, &v)| v).count() as f64;

    m1 * q.ln() + m0 * (1. - q).ln()
}

pub fn log_likelihood(
    alpha: f64,
    beta: f64,
    obs: &HashMap<GeneSymbol, bool>,
    hid: &HashMap<GeneSymbol, bool>,
) -> f64 {
    let [n00, n01, n10, n11] = counts_nxy(obs, hid);
    (n00 as f64) * (1.0 - alpha).ln()
        + (n10 as f64) * alpha.ln()
        + (n01 as f64) * beta.ln()
        + (n11 as f64) * (1.0 - beta).ln()
}

pub fn log_posterior_unnorm(
    alpha: f64,
    beta: f64,
    q: f64,
    obs: &HashMap<GeneSymbol, bool>,
    hid: &HashMap<GeneSymbol, bool>,
    terms: &HashMap<TermId, bool>,
) -> f64 {
    log_likelihood(alpha, beta, obs, hid) + log_prior(q, terms)
}
