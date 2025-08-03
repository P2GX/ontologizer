use crate::calculation::results::AnalysisResults;

pub trait MultipleTestingCorrection {
    fn adjust_pvalues(&self, results: &mut AnalysisResults);
    fn name(&self) -> &'static str;
}
