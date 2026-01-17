use std::collections::HashSet;
use ontolius::ontology::csr::FullCsrOntology;
use crate::core::AnnotationIndex;
use crate::core::result::EnrichmentResult;

pub fn run(ontology : &FullCsrOntology, annotation_index: AnnotationIndex, study_genes : HashSet<String>) -> EnrichmentResult {
    let n_population_genes = annotation_index.get_genes().len();
    let n_study_genes = study_genes.len();

    let genes_to_terms = annotation_index.get_genes_to_terms(true);

    todo!();
    // let study_genes : Vec<bool> = (0..n_genes)
    //     .map(|i| study_genes.contains(annotation_index.get_gene_by_index(i)))
    //     .collect();
    //
    // let test = Hypergeometric::new();
    //
    // let measure : PValue = test.calculate();
    //
    // let mut result = EnrichmentResult::from_measure(
    //     &measure,
    //     &ontology,
    //     annotation_index.get_terms(),
    //     annotation_index.get_genes(),
    //     &observed_genes,
    //     &terms_to_genes,
    // );
    //
    // // Optional: Sort
    // result.sort_by_score(true); // descending for probability
    //
    // result
}
