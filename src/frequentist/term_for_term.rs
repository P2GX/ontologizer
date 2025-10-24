use std::collections::HashSet;

use crate::frequentist::hypergeo::Hypergeometric;

use ontolius::{
    ontology::{csr::FullCsrOntology, OntologyTerms},
    term::MinimalTerm,
};
use crate::core::AnnotationIndex;
use crate::core::{GeneSet, GeneSymbol};
use crate::frequentist::PValueCalculation;
use super::results::{get_term_aspect, AnalysisResults, GOTermResult};

pub struct TermForTerm;
#[allow(non_snake_case)]
impl PValueCalculation for TermForTerm {
    fn calculate_p_values(
        &self,
        go: &FullCsrOntology,
        annotation_container: &AnnotationIndex,
        study: &GeneSet,
        population: &GeneSet,
        results: &mut AnalysisResults,
    ) -> () {
        let study_genes = study.recognized_genes();
        let population_genes = population.recognized_genes();
        let n = study_genes.len();
        let N = population_genes.len();

        let study_terms_count = annotation_container.term_counts_for_subset(study, go);
        let population_terms_count = annotation_container.term_counts_for_subset(population, go);

        for (term, &k) in study_terms_count.iter() {
                let &K = population_terms_count.get(term).unwrap();

                // let K = annotated_genes.len();
                let mut hypergeom = Hypergeometric::new();

                let raw_p_value;
                if k > 1 {
                    raw_p_value = hypergeom.phyper(k - 1, n, K, N, false).unwrap();
                }
                else {
                    continue;
                    // todo!(this was previously skipped, correctly the p-value should be 1.0)
                    raw_p_value = 1.0;
                }

                let term_id = go.term_by_id(term);
                let result = GOTermResult::new(
                    term_id.unwrap().name().to_string(),
                    term.to_string(),
                    (k as u32, K as u32, n as u32, N as u32),
                    raw_p_value as f32,
                    raw_p_value as f32, // use raw pvalue as default
                    get_term_aspect(go, term),
                );

                results.add_result(result);
            }
        results.sort_by_p_value();
    }
}
// Computes the intersection of a gene set with the set of genes annotated to a specific GO term.
fn gene_set_intersection(
    gene_set: &HashSet<GeneSymbol>,
    genes_annotated_to_term: &HashSet<GeneSymbol>,
) -> HashSet<GeneSymbol> {
    let intersection: HashSet<_> = gene_set
        .intersection(genes_annotated_to_term)
        .cloned()
        .collect();
    intersection
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::frequentist::mtc::{Bonferroni, MultipleTestingCorrection};

    use crate::frequentist::results::{AnalysisResults, MethodEnum, MtcEnum};
    use crate::core::{load_gene_set, separate_gene_set};
    use crate::core::Ontologizer;

    #[test]
    // #[ignore]
    fn test_go_analysis() {
        // Paths to the necessary files
        // These paths should be adjusted based on your local setup
        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let pop_set_path = "tests/data/population.txt";
        let study_set_path = "tests/data/study.txt";

        let mtc_method = MtcEnum::Bonferroni;
        // Load the GO ontology
        let go = Ontologizer::new(go_path);
        let go_ref = go.ontology();

        // Load the GOA annotations
        let mut annotation_container = AnnotationIndex::new(gaf_path, go_ref);

        // Load the population and study gene sets
        let study_gene_symbols = load_gene_set(study_set_path).expect("Failed to parse study gene set");
        let study_gene_set = separate_gene_set(&annotation_container.get_annotations(), study_gene_symbols);

        let pop_gene_symbols = load_gene_set(pop_set_path).expect("Failed to parse population gene set");
        let pop_gene_set = separate_gene_set(&annotation_container.get_annotations(), pop_gene_symbols);

        let mut results = AnalysisResults::new(MethodEnum::TermForTerm, mtc_method);

        TermForTerm.calculate_p_values(go_ref, &annotation_container, &study_gene_set, &pop_gene_set, &mut results);
        Bonferroni.adjust_pvalues(&mut results);
        results
            .write_tsv("tests/data/enrichment_results.tsv")
            .expect("Failed to write results");
    }
}
