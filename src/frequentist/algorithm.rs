use crate::frequentist::hypergeometric::Hypergeometric;
use indexmap::IndexMap;

use super::results::{get_term_aspect, AnalysisResults, GOTermResult};
use crate::core::AnnotationIndex;
use crate::core::GeneSet;
use crate::frequentist::PValueCalculation;
use ontolius::{
    ontology::{csr::FullCsrOntology, OntologyTerms},
    term::MinimalTerm,
    TermId,
};

pub struct TermForTerm;
#[allow(non_snake_case)]
impl PValueCalculation for TermForTerm {
    fn calculate_p_values(
        &self,
        go: &FullCsrOntology,
        annotation_index: &AnnotationIndex,
        study: &GeneSet,
        population: &GeneSet,
        results: &mut AnalysisResults,
    ) -> () {
        let study_genes = study.recognized_genes();
        let population_genes = population.recognized_genes();
        let n = study_genes.len();
        let N = population_genes.len();

        let study_terms_count = term_count_for_subset(study, annotation_index);
        let population_terms_count = term_count_for_subset(population, annotation_index);

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
    let mut dense_counts = vec![0usize; annotations.get_terms().len()];

    for gene in gene_set.recognized_genes() {
        if let Some(gene_idx) = annotations.get_index_by_gene(gene) {
            for &term_idx in annotations.get_term_idxs_for_gene_idx(gene_idx) {
                dense_counts[term_idx] += 1;
            }
        }
    }

    let mut counts = IndexMap::new();
    for (term_idx, &count) in dense_counts.iter().enumerate() {
        if count > 0 {
            // Safe unwrap: term_map and dense_counts are perfectly aligned by design
            let term_id = annotations.get_term_by_index(term_idx);
            counts.insert(term_id.clone(), count);
        }
    }

    counts
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::core::{load_gene_set, separate_gene_set};
    use crate::frequentist::correction::Bonferroni;
    use crate::frequentist::results::{AnalysisResults, MethodEnum, MtcEnum};
    use oboannotation::go::stats::get_annotation_map;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use std::process;
    use crate::frequentist::correction::MultipleTestingCorrection;

    #[test]
    fn test_go_analysis() {
        // Paths to the necessary files
        // These paths should be adjusted based on your local setup
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let population_genes_path = "tests/data/GOnone/population.txt";
        let study_genes_path = "tests/data/GOnone/study.txt";

        let mtc_method = MtcEnum::Bonferroni;

        // ------ Load Gene Sets ------
        let study_genes = load_gene_set(study_genes_path).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        let population_genes = load_gene_set(population_genes_path).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        // ------ Load Gene Ontology and Annotations ------
        let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
        let ontology: FullCsrOntology =
            ontology_loader
                .load_from_path(go_path)
                .unwrap_or_else(|err| {
                    eprintln!("Error: {}", err);
                    process::exit(1);
                });

        let annotations_loader = GoGafAnnotationLoader;
        let annotations: GoAnnotations = annotations_loader
            .load_from_path(gaf_path)
            .expect("Could not load GAF file");

        let annotation_index =
            AnnotationIndex::new(annotations, &ontology, Some(&population_genes));
        let annotated_genes = get_annotation_map(&annotation_index.annotations)
            .into_keys()
            .collect();

        let study_gene_set = separate_gene_set(&annotated_genes, study_genes);
        let pop_gene_set = separate_gene_set(&annotated_genes, population_genes);

        let mut results = AnalysisResults::new(MethodEnum::TermForTerm, mtc_method);

        TermForTerm.calculate_p_values(
            &ontology,
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
