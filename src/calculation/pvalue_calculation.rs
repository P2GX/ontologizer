use crate::annotations::AnnotationContainer;
use crate::calculation::results::AnalysisResults;
use crate::gene_set::GeneSet;
use ontolius::ontology::csr::FullCsrOntology;
pub trait PValueCalculation {
    fn calculate_p_values(
        &self,
        go: &FullCsrOntology,
        annotation_container: &AnnotationContainer,
        study: &GeneSet,
        population: &GeneSet,
         results: &mut AnalysisResults,
    ) -> ();
}
