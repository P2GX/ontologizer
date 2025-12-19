#[cfg(test)]
mod test {
    use crate::bayesian::algorithm::{Algorithm, MetropolisHasting};
    use crate::bayesian::model::OrModel;
    use crate::bayesian::proposer::UniformProposer;
    use crate::bayesian::recorder::Frequency;

    use crate::bayesian::state::MgsaState;
    use crate::core::{AnnotationIndex, Ontologizer, load_gene_set};

    use crate::core::result::{AnalysisResult, BayesianResult};
    use csv::Writer;

    #[test]
    fn test_mgsa() {
        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let study_set_path = "tests/data/study.txt";

        // Load the GO ontology
        let go = Ontologizer::new(go_path);
        let go_ref = go.ontology();
        let annotations = AnnotationIndex::new(gaf_path, go_ref);
        // Load the GOA annotations

        // Load the population and study gene sets
        let gene_symbols = load_gene_set(study_set_path).expect("Failed to parse study gene set");

        // MGSA Parameter
        let p = 0.5;
        let alpha = 0.05;
        let beta = 0.10;

        let terms_to_genes = annotations.get_terms_to_genes(true);
        let n_genes = annotations.genes().len();
        let observed_genes: Vec<bool> = (0..n_genes)
            .map(|i| gene_symbols.contains(annotations.get_index_gene(i)))
            .collect();

        let model = OrModel::new(terms_to_genes.clone(), observed_genes, beta, p, alpha);
        let terms = model.heuristic_start();

        let mut state = MgsaState::new(terms);

        let proposer = UniformProposer::new();

        let mut algorithm = MetropolisHasting::new(model, proposer, 1_000_000, 100_000);

        let measure: Frequency = algorithm.sample(&mut state);

        // Create the Result (Eagerly resolves all strings)
        let mut result = BayesianResult::from_counts(
            &measure,
            &go_ref,
            annotations.terms(),
            annotations.genes(),
            &terms_to_genes,
            0.5,
            0.05,
            0.05,
        );

        // Optional: Sort
        result.sort_by_score(true); // descending for probability

        // Serialize to CSV
        let mut wtr = Writer::from_path("results.csv").unwrap();
        for item in result.items() {
            wtr.serialize(item).unwrap();
        }
    }
}
