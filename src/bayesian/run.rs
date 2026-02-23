use crate::bayesian::algorithm::{Algorithm, MetropolisHastings};
use crate::bayesian::model::OrModel;
use crate::bayesian::proposer::{MixedProposer, ParameterGaussProposer, TermToggleSwapProposer};
use crate::bayesian::recorder::MgsaRecorder;
use crate::bayesian::state::MgsaState;
use crate::core::AnnotationIndex;
use crate::core::result::EnrichmentResult;
use ontolius::ontology::csr::FullCsrOntology;
use std::collections::HashSet;

pub fn analysis(
    ontology: &FullCsrOntology,
    annotations: AnnotationIndex,
    study_genes: HashSet<String>,
) -> EnrichmentResult {
    // --- Configuration ---
    let p = 0.05;
    let alpha = 0.05;
    let beta = 0.10;

    let iterations = 50_000_000;
    let burn_in = 1_000_000;
    let term_update_prob = 0.95; // 95% Term moves, 5% Parameter moves

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
    let model = OrModel::new(terms_to_genes.clone(), obs_genes.clone(), p, alpha, beta);
    let mut state = MgsaState::new(vec![false; n_terms], p, alpha, beta);

    // --- Algorithm Setup ---
    let proposer = MixedProposer::new(
        TermToggleSwapProposer::new(),
        ParameterGaussProposer::new(1.0),
        term_update_prob,
    );
    let mut algorithm = MetropolisHastings::new(model, proposer, iterations, burn_in);

    // --- Execution ---
    let (measures, parameter) = algorithm.sample::<MgsaRecorder>(&mut state);

    // --- Result Generation ---
    let mut result =
        EnrichmentResult::from_measures(&measures, &ontology, &annotations, &obs_genes);

    result.sort_by_score(true); // Sort descending by probability

    result
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

        (ontology, annotation_index, study_genes.recognized_genes())
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

        let result = super::analysis(&ontology, annotation_index, study_genes);

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

        let result = super::analysis(&ontology, annotation_index, study_genes);

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
