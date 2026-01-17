mod algorithm;
mod correction;
pub mod hypergeometric;
mod measure;

mod results;
mod run;

pub use run::run;

use ontolius::ontology::csr::FullCsrOntology;
pub use correction::{MultipleTestingCorrection};
use crate::core::{AnnotationIndex, GeneSet};
use crate::frequentist::results::AnalysisResults;

pub trait PValueCalculation {
    fn calculate_p_values(
        &self,
        go: &FullCsrOntology,
        annotation_container: &AnnotationIndex,
        study: &GeneSet,
        population: &GeneSet,
        results: &mut AnalysisResults,
    ) -> ();
}

