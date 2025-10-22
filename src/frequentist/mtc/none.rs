use super::MultipleTestingCorrection;
use crate::frequentist::results::AnalysisResults;

pub struct None;

impl MultipleTestingCorrection for None {
    fn adjust_pvalues(&self, _analysis_results: &mut AnalysisResults) {
        // No adjustment is made
    }

    fn name(&self) -> &'static str {
        "No Correction"
    }
}
