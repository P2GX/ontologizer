use crate::frequentist::distribution::DiscreteDistribution;

pub trait StatisticalTest {
    /// The input type for the test (e.g., the observed count).
    type Input;

    /// Calculates the p-value given a distribution and an observation.
    fn calculate<D: DiscreteDistribution>(&self, distribution: &D, observed: Self::Input) -> f64;
}

/// A standard one-sided test (enrichment). Checks if we have more overlap than expected: P(X >= k)
pub struct OneSidedEnrichmentTest;

impl StatisticalTest for OneSidedEnrichmentTest {
    type Input = usize;

    fn calculate<D: DiscreteDistribution>(&self, dist: &D, k: usize) -> f64 {
        // P(X >= k) is equivalent to P(X > k - 1)
        if k == 0 {
            1.0
        } else {
            dist.sf(k - 1)
        }
    }
}