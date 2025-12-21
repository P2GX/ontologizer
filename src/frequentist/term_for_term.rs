use crate::frequentist::hypergeo::Hypergeometric;
use indexmap::IndexMap;

use super::results::{AnalysisResults, GOTermResult, get_term_aspect};
use crate::core::AnnotationIndex;
use crate::core::GeneSet;
use crate::frequentist::PValueCalculation;
use ontolius::{
    TermId,
    ontology::{OntologyTerms, csr::FullCsrOntology},
    term::MinimalTerm,
};

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

        let study_terms_count = term_count_for_subset(study, annotation_container);
        let population_terms_count = term_count_for_subset(population, annotation_container);

        for (term, &k) in study_terms_count.iter() {
            // TODO
            let &K = population_terms_count.get(term).unwrap();

            // let K = annotated_genes.len();
            let mut hypergeom = Hypergeometric::new();

            let raw_p_value;
            if k > 1 {
                raw_p_value = hypergeom.phyper(k - 1, n, K, N, false).unwrap();
            } else {
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

pub fn term_count_for_subset(
    gene_set: &GeneSet,
    annotations: &AnnotationIndex,
) -> IndexMap<TermId, usize> {
    let mut dense_counts = vec![0usize; annotations.terms().len()];

    for gene in gene_set.recognized_genes() {
        if let Some(gene_idx) = annotations.get_gene_index(gene) {
            for &term_idx in annotations.get_term_idxs_for_gene_idx(gene_idx) {
                dense_counts[term_idx] += 1;
            }
        }
    }

    let mut counts = IndexMap::new();
    for (term_idx, &count) in dense_counts.iter().enumerate() {
        if count > 0 {
            // Safe unwrap: term_map and dense_counts are perfectly aligned by design
            let term_id = annotations.get_index_term(term_idx);
            counts.insert(term_id.clone(), count);
        }
    }

    counts
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::frequentist::mtc::{Bonferroni, MultipleTestingCorrection};
    use oboannotation::go::stats::get_annotation_map;

    use crate::core::Ontologizer;
    use crate::core::{load_gene_set, separate_gene_set};
    use crate::frequentist::results::{AnalysisResults, MethodEnum, MtcEnum};

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

        // Load the population and study gene sets
        let study_gene_symbols =
            load_gene_set(study_set_path).expect("Failed to parse study gene set");

        let pop_gene_symbols =
            load_gene_set(pop_set_path).expect("Failed to parse population gene set");

        // Load the GOA annotations
        let annotation_index = AnnotationIndex::new(gaf_path, go_ref, &pop_gene_symbols);
        let annotated_genes = get_annotation_map(&annotation_index.annotations)
            .into_keys()
            .collect();

        let study_gene_set = separate_gene_set(&annotated_genes, study_gene_symbols);
        let pop_gene_set = separate_gene_set(&annotated_genes, pop_gene_symbols);

        let mut results = AnalysisResults::new(MethodEnum::TermForTerm, mtc_method);

        TermForTerm.calculate_p_values(
            go_ref,
            &annotation_index,
            &study_gene_set,
            &pop_gene_set,
            &mut results,
        );
        Bonferroni.adjust_pvalues(&mut results);
        results
            .write_tsv("tests/data/enrichment_results.tsv")
            .expect("Failed to write results");
    }
}
