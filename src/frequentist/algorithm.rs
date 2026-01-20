use crate::frequentist::distribution::DiscreteDistribution;

pub trait StatisticalTest<D>
where
    D: DiscreteDistribution,
{
    // Returns the raw p-value
    fn calculate(&self, distribution: &D, observed: usize) -> f64;
}

pub struct OneSidedEnrichmentTest;

impl<D> StatisticalTest<D> for OneSidedEnrichmentTest
where
    D: DiscreteDistribution,
{
    fn calculate(&self, dist: &D, k: usize) -> f64 {
        if k == 0 { 1.0 } else { dist.sf(k - 1) }
    }
}
