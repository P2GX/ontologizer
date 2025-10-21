use crate::annotations::AnnotationIndex;
use crate::calculation::results::AnalysisResults;
use crate::geneset::GeneSet;
use ontolius::ontology::csr::FullCsrOntology;
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
