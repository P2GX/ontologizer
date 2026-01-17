use crate::frequentist::results::AnalysisResults;

pub trait MultipleTestingCorrection {
    fn adjust_pvalues(&self, results: &mut AnalysisResults);
    fn name(&self) -> &'static str;
}

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

pub struct BonferroniHolm;

impl MultipleTestingCorrection for BonferroniHolm {
    fn adjust_pvalues(&self, analysis_results: &mut AnalysisResults) {
        analysis_results.sort_by_p_value();
        let n: f32 = analysis_results.num_hypotheses();
        let mut i: f32 = 1.0;
        for result in analysis_results.iter_mut() {
            result.set_adj_pval(result.p_val() * (n - i + 1.0));
            i += 1.0;
        }

        enforce_pvalue_monotony(analysis_results);
    }

    fn name(&self) -> &'static str {
        "Bonferroni-Holm"
    }
}

fn enforce_pvalue_monotony(analysis_results: &mut AnalysisResults) {
    // let n: f32 = analysis_results.num_hypotheses();
    let mut prev: f32 = 1.0;
    for result in analysis_results.iter_mut().rev() {
        let curr = result.adj_pval();
        let new_val = curr.min(prev);
        result.set_adj_pval(new_val);
        prev = new_val;
    }
}

pub struct None;

impl MultipleTestingCorrection for None {
    fn adjust_pvalues(&self, _analysis_results: &mut AnalysisResults) {
        // No adjustment is made
    }

    fn name(&self) -> &'static str {
        "No Correction"
    }
}