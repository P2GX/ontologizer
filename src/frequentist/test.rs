//
// use crate::frequentist::distribution::{DiscreteDistribution};
//
// pub trait StatisticalTest<D>
// where
//     D: DiscreteDistribution,
// {
//     /// Calculates the measure (e.g., P-value) given a specific distribution instance
//     /// and the observed data.
//     fn calculate(&self, distribution: &D, observed: usize) -> PValue;
// }
//
// /// A standard one-sided test (enrichment). Checks if we have more overlap than expected: P(X >= k)
// pub struct OneSidedEnrichmentTest;
//
// impl<D> StatisticalTest<D> for OneSidedEnrichmentTest
// where
//     D: DiscreteDistribution,
// {
//     fn calculate(&self, dist: &D, k: usize) -> PValue {
//         // P(X >= k) = P(X > k - 1) = sf(k - 1)
//         let p_val = if k == 0 { 1.0 } else { dist.sf(k - 1) };
//
//         PValue::from_single(p_val)
//     }
// }
