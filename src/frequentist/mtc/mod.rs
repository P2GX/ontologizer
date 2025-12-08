use crate::frequentist::results::AnalysisResults;

pub trait MultipleTestingCorrection {
    fn adjust_pvalues(&self, results: &mut AnalysisResults);
    fn name(&self) -> &'static str;
}

mod bonferroni;
mod bonferroni_holm;
mod none;

pub use bonferroni::Bonferroni;
pub use bonferroni_holm::BonferroniHolm;
pub use none::None;
