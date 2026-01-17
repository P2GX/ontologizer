use std::collections::HashSet;
use ontolius::ontology::csr::FullCsrOntology;
use crate::bayesian::algorithm::{Algorithm, MetropolisHasting};
use crate::bayesian::measure::Probability;
use crate::bayesian::model::OrModel;
use crate::bayesian::proposer::UniformToggleProposer;
use crate::bayesian::recorder::ProbabilityRecorder;
use crate::bayesian::state::MgsaState;
use crate::core::AnnotationIndex;
use crate::core::result::EnrichmentResult;

pub fn run(ontology : &FullCsrOntology, annotation_index: AnnotationIndex, study_genes : HashSet<String>) -> EnrichmentResult {
    let p = 0.01;
    let alpha = 0.05;
    let beta = 0.10;

    let n_genes = annotation_index.get_genes().len();
    let n_terms = annotation_index.get_terms().len();

    let terms_to_genes = annotation_index.get_terms_to_genes(true);

    let observed_genes : Vec<bool> = (0..n_genes)
        .map(|i| study_genes.contains(annotation_index.get_gene_by_index(i)))
        .collect();

    let model = OrModel::new(
        terms_to_genes.clone(),
        observed_genes.clone(),
        p,
        alpha,
        beta,
    );

    let mut state = MgsaState::new(vec![false; n_terms]);
    let proposer = UniformToggleProposer::new();
    let mut algorithm = MetropolisHasting::new(model, proposer, 50_000_000, 1_000_000);

    let measure: Probability = algorithm.sample::<ProbabilityRecorder>(&mut state);

    // Create the Result (Eagerly resolves all strings)
    let mut result = EnrichmentResult::from_measure(
        &measure,
        &ontology,
        annotation_index.get_terms(),
        annotation_index.get_genes(),
        &observed_genes,
        &terms_to_genes,
    );

    // Optional: Sort
    result.sort_by_score(true); // descending for probability
    
    result
}


#[cfg(test)]
mod test {
    use super::*; // Use super to access run()
    use crate::core::load_gene_set;
    use crate::core::AnnotationIndex;
    use csv::Writer;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use ontolius::ontology::csr::FullCsrOntology;
    use std::process;
    #[test]
    fn test_mgsa() {
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let study_set_path = "tests/data/GOnone/study.txt";
        let pop_set_path = "tests/data/GOnone/population.txt";

        // Load the population and study gene sets
        let study_genes = load_gene_set(study_set_path).expect("Failed to parse study genes");
        let population_genes =
            load_gene_set(pop_set_path).expect("Failed to parse population genes");

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

        let result = super::run(&ontology, annotation_index, study_genes);

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GOnone/results.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }

    #[test]
    fn test_specific_term() {
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let study_set_path = "tests/data/GO0090717/study.txt";
        let pop_set_path = "tests/data/GO0090717/population.txt";

        // Load the population and study gene sets
        let study_genes = load_gene_set(study_set_path).expect("Failed to parse study gene set");
        let population_genes = load_gene_set(pop_set_path).expect("Failed to parse study gene set");

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

        let result = super::run(&ontology, annotation_index, study_genes);

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GO0090717/results.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }
}
