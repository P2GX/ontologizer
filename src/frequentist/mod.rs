use crate::core::{AnnotationIndex, GeneSet};
use crate::frequentist::results::AnalysisResults;
use ontolius::ontology::csr::FullCsrOntology;

mod hypergeo;
pub mod mtc;
pub mod results; // Contains the results of the enrichment analysis
pub mod term_for_term;
// Contains the term-for-term enrichment analysis

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
