use ontolius::ontology::csr::FullCsrOntology;
use crate::core::{AnnotationIndex, GeneSet};
use crate::frequentist::results::AnalysisResults;

pub mod results; // Contains the results of the enrichment analysis
pub mod term_for_term;
pub mod mtc;
mod hypergeom;
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