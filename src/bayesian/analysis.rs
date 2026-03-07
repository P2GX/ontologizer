use crate::bayesian::algorithm::{Algorithm, MetropolisHastings};
use crate::bayesian::model::OrModel;
use crate::bayesian::proposer::{MixedProposer, ParameterGaussProposer, TermToggleSwapProposer};
use crate::bayesian::recorder::MgsaRecorder;
use crate::bayesian::state::MgsaState;
use crate::core::AnnotationIndex;
use crate::core::result::{AnalysisResult, Measure};
use indexmap::IndexSet;
use ontolius::ontology::csr::FullCsrOntology;
use std::collections::HashSet;

pub fn approximate_minimal_term_cover(
    gene_idxs: &IndexSet<usize>,
    term_idxs: &IndexSet<usize>,
    annotations: &AnnotationIndex,
) -> IndexSet<usize> {
    let mut uncovered_genes: IndexSet<usize> = gene_idxs.clone();
    let mut available_terms: IndexSet<usize> = term_idxs.clone();
    let mut minimal_terms: IndexSet<usize> = IndexSet::new();

    while !uncovered_genes.is_empty() {
        let mut best_term: Option<usize> = None;
        let mut max_covered: usize = 0;

        for &term_idx in &available_terms {
            let term_genes: &IndexSet<usize> = annotations.get_gene_idxs_for_term_idx(term_idx);
            let covered_count: usize = term_genes.intersection(&uncovered_genes).count();

            if covered_count > max_covered {
                max_covered = covered_count;
                best_term = Some(term_idx);
            }
        }

        match best_term {
            Some(term_idx) => {
                minimal_terms.insert(term_idx);
                available_terms.swap_remove(&term_idx);

                let term_genes: &IndexSet<usize> = annotations.get_gene_idxs_for_term_idx(term_idx);
                for &g in term_genes {
                    // swap_remove provides O(1) removal by swapping the removed element
                    // with the last element, sacrificing order for performance.
                    uncovered_genes.swap_remove(&g);
                }
            }
            None => break,
        }
    }

    minimal_terms
}

pub fn analysis(
    ontology: &FullCsrOntology,
    annotations: &AnnotationIndex,
    study_genes: &HashSet<String>,
) -> AnalysisResult {
    // --- Configuration ---

    // Init prior probability for term activation by assessing the minimal set of terms required
    // to explain observed genes
    let n_terms = annotations.get_terms().len();
    let gene_idxs: IndexSet<usize> = study_genes
        .iter()
        .filter_map(|g| annotations.get_idx_by_gene(g))
        .collect();

    let term_idxs: IndexSet<usize> = gene_idxs
        .iter()
        .flat_map(|&g| annotations.get_term_idxs_for_gene_idx(g).iter().copied())
        .collect();
    let n_min_terms = approximate_minimal_term_cover(&gene_idxs, &term_idxs, annotations).len();

    let p_init = n_min_terms as f64 / n_terms as f64;
    let alpha_init = 0.05;
    let beta_init = 0.05;

    let iterations = annotations.get_terms().len() * 5_000;
    let burn_in = annotations.get_terms().len() * 1_000;

    // --- Data Preparation ---
    let n_genes = annotations.get_genes().len();
    let n_terms = annotations.get_terms().len();

    let terms_to_genes = annotations.get_terms_to_genes(true);

    // Map string-based study genes to the internal boolean observation vector.
    // vec[i] = true if the gene at index i is in the study set.
    let obs_genes: Vec<bool> = (0..n_genes)
        .map(|i| study_genes.contains(annotations.get_gene_by_idx(i)))
        .collect();

    // --- Model & State Initialization ---
    let model = OrModel::new(
        terms_to_genes.clone(),
        obs_genes.clone(),
        p_init,
        alpha_init,
        beta_init,
    );
    let mut state = MgsaState::new(vec![false; n_terms], p_init, alpha_init, beta_init);

    // --- Algorithm Setup ---
    let proposer = MixedProposer::new(
        TermToggleSwapProposer::new(),
        ParameterGaussProposer::new(1.0),
        0.90, // 90% Term moves, 10% Parameter moves
    );
    let mut algorithm = MetropolisHastings::new(model, proposer, iterations, burn_in);

    // --- Execution ---
    let (measures, parameter) = algorithm.sample::<MgsaRecorder>(&mut state);

    // --- Result Generation ---
    let mut result = AnalysisResult::from_measures(&measures, &ontology, &annotations, &obs_genes);

    result.sort_by_score(true); // Sort descending by probability

    let p_str = format!("{:.4}", parameter[0].score());
    let alpha_str = format!("{:.4}", parameter[1].score());
    let beta_str = format!("{:.4}", parameter[2].score());

    result.with_meta(&[
        ("Method", "Bayesian"),
        ("p", &p_str),
        ("alpha", &alpha_str),
        ("beta", &beta_str),
    ])
}

#[cfg(test)]
mod test {
    use crate::GeneSet;
    use crate::core::AnnotationIndex;
    use csv::Writer;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use ontolius::ontology::csr::FullCsrOntology;
    use std::collections::HashSet;
    use std::path::Path;
    use std::process;

    /// Helper to load test data and avoid code duplication in tests.
    fn load_test_data(
        go_path: &str,
        gaf_path: &str,
        study_path: &str,
        pop_path: &str,
    ) -> (FullCsrOntology, AnnotationIndex, HashSet<String>) {
        // Verify files exist before attempting to load (avoids confusing panic messages)
        if !Path::new(go_path).exists() {
            panic!(
                "Test file missing: {}. Please ensure test data is present.",
                go_path
            );
        }

        let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
        let ontology: FullCsrOntology = ontology_loader
            .load_from_path(go_path)
            .expect("Failed to load Ontology");

        let annotations_loader = GoGafAnnotationLoader;
        let annotations: GoAnnotations = annotations_loader
            .load_from_path(gaf_path)
            .expect("Failed to load GAF");

        let study_genes = GeneSet::from_file(study_path, &annotations).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        let population_genes = GeneSet::from_file(pop_path, &annotations).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        // Construct the index
        let annotation_index = AnnotationIndex::new(
            annotations,
            &ontology,
            Some(&population_genes.recognized_genes()),
        );

        (
            ontology,
            annotation_index,
            study_genes.recognized_genes().clone(),
        )
    }

    #[test]
    #[ignore = "Long running integration test requires external data"]
    fn test_mgsa_execution() {
        let (ontology, annotation_index, study_genes) = load_test_data(
            "tests/data/GO/go-basic.json",
            "tests/data/GO/goa_human.gaf",
            "tests/data/GOnone/study.txt",
            "tests/data/GOnone/population.txt",
        );

        let result = super::analysis(&ontology, &annotation_index, &study_genes);

        // Assertions: Verify we got results and they are valid probabilities
        assert!(!result.items.is_empty(), "Result set should not be empty");
        assert!(result.items[0].score <= 1.0 && result.items[0].score >= 0.0);
    }

    #[test]
    #[ignore = "Long running integration test requires external data"]
    fn test_specific_term_recovery() {
        // This test uses a study set generated specifically from term GO:0090717.
        // We expect this term to be discovered with high probability.
        let target_term = "GO:0090717";

        let (ontology, annotation_index, study_genes) = load_test_data(
            "tests/data/GO/go-basic.json",
            "tests/data/GO/goa_human.gaf",
            "tests/data/GO0090717/study.txt",
            "tests/data/GO0090717/population.txt",
        );

        let result = super::analysis(&ontology, &annotation_index, &study_genes);

        // Check if the target term is in the top results
        let target_item = result.items.iter().find(|i| i.id == target_term);

        match target_item {
            Some(item) => {
                println!(
                    "Found target {} with probability {:.4}",
                    target_term, item.score
                );
            }
            None => {
                panic!("Target term {} was not found in the results", target_term);
            }
        }

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GO0090717/results.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }
}
