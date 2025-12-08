use crate::frequentist::results::AnalysisResults;

pub trait MultipleTestingCorrection {
    fn adjust_pvalues(&self, results: &mut AnalysisResults);
    fn name(&self) -> &'static str;
}

mod none;
mod  bonferroni;
mod bonferroni_holm;

pub use none::None;
pub use bonferroni::Bonferroni;
pub use bonferroni_holm::BonferroniHolm;



