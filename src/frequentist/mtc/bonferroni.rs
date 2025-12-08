use super::MultipleTestingCorrection;
use crate::frequentist::results::AnalysisResults;

pub struct Bonferroni;

impl MultipleTestingCorrection for Bonferroni {
    fn adjust_pvalues(&self, analysis_results: &mut AnalysisResults) {
        let n = analysis_results.num_hypotheses();
        for result in analysis_results.iter_mut() {
            result.set_adj_pval((result.p_val() * n).min(1.0));
        }
    }

    fn name(&self) -> &'static str {
        "Bonferroni"
    }
}
