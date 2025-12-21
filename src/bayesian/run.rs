#[cfg(test)]
mod test {
    use crate::bayesian::algorithm::{Algorithm, MetropolisHasting};
    use crate::bayesian::model::OrModel;
    use crate::bayesian::proposer::UniformProposer;
    use crate::bayesian::recorder::Frequency;

    use crate::bayesian::state::MgsaState;
    use crate::core::{AnnotationIndex, Ontologizer, load_gene_set, separate_gene_set};

    use crate::core::result::{AnalysisResult, BayesianResult};
    use csv::Writer;

    use std::time::Instant;

    fn approx_equal(a: f64, b : f64, epsilon : f64) -> bool {
        (a-b).abs() < epsilon
    }

    #[test]
    fn test_mgsa() {
        let start_total = Instant::now();
        let mut last_checkpoint = Instant::now();

        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let study_set_path = "tests/data/study.txt";
        let pop_set_path = "tests/data/population.txt";

        // Load the population and study gene sets
        let obs_gene_symbols =
            load_gene_set(study_set_path).expect("Failed to parse study gene set");
        let pop_gene_symbols =
            load_gene_set(pop_set_path).expect("Failed to parse population gene set");
        println!("Gene Set Load: {:.2?}", last_checkpoint.elapsed());
        last_checkpoint = Instant::now(); // Reset timer

        // Load the GO ontology
        let go = Ontologizer::new(go_path);
        let go_ref = go.ontology();
        println!("Ontology Load: {:.2?}", last_checkpoint.elapsed());
        last_checkpoint = Instant::now();

        // Build the AnnotationIndex restricted to Population Gene Set.
        let annotations = AnnotationIndex::new(gaf_path, go_ref, &pop_gene_symbols);
        println!("Annotations Load: {:.2?}", last_checkpoint.elapsed());
        last_checkpoint = Instant::now();

        // MGSA Parameter
        let p = 0.5;
        let alpha = 0.05;
        let beta = 0.10;

        let terms_to_genes = annotations.get_terms_to_genes(true);
        let n_genes = annotations.genes().len();
        let observed_genes: Vec<bool> = (0..n_genes)
            .map(|i| obs_gene_symbols.contains(annotations.get_index_gene(i)))
            .collect();

        let model = OrModel::new(
            terms_to_genes.clone(),
            observed_genes.clone(),
            p,
            alpha,
            beta,
        );
        let terms = model.heuristic_start();
        let mut state = MgsaState::new(terms);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 5_000_000, 1_000_000);

        println!("Model Setup: {:.2?}", last_checkpoint.elapsed());
        last_checkpoint = Instant::now();

        let measure: Frequency = algorithm.sample(&mut state);
        println!("MCMC Sampling: {:.2?}", last_checkpoint.elapsed());
        last_checkpoint = Instant::now();

        // Create the Result (Eagerly resolves all strings)
        let mut result = BayesianResult::from_counts(
            &measure,
            &go_ref,
            annotations.terms(),
            annotations.genes(),
            &observed_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
        );

        // Optional: Sort
        result.sort_by_score(true); // descending for probability

        // Serialize to CSV
        let mut wtr = Writer::from_path("results.csv").unwrap();
        for item in result.items() {
            wtr.serialize(item).unwrap();
        }

        println!(
            "Result Processing & CSV Write: {:.2?}",
            last_checkpoint.elapsed()
        );

        println!(
            "Result Processing & CSV Write: {:.2?}",
            last_checkpoint.elapsed()
        );
        println!("Total Execution Time: {:.2?}", start_total.elapsed());
    }

    #[test]
    fn test_one_term_one_gene(){
        let state = vec![false];
        let observations = vec![true];
        let state_to_observations : Vec<Vec<usize>> = vec![
            vec![0]
        ];

        let p = 0.5;
        let alpha = 0.05;
        let beta = 0.1;

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Frequency = algorithm.sample(&mut state);
        println!{"{:?}", measure.frequencies}

        let post_sim = measure.frequencies[0];
        let post_theory = p*(1.-beta)/(alpha*(1.-p) + (1.-beta)*p);
        assert!(approx_equal(post_sim, post_theory, 0.01));
    }

    #[test]
    fn test_two_term_one_gene(){
        let state = vec![false, false];
        let observations = vec![true];
        let state_to_observations : Vec<Vec<usize>> = vec![
            vec![0],
            vec![0]
        ];

        let p = 0.5;
        let beta = 0.1;
        let alpha = 0.05;

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Frequency = algorithm.sample(&mut state);
        println!{"{:?}", measure.frequencies}

        let post_t1_sim = measure.frequencies[0];
        let post_t2_sim = measure.frequencies[1];
        let post_theory = p*(1.-beta)/(alpha*(1.-p).powf(2.) + (1.-beta)*(1.-(1.-p).powf(2.)));
        assert!(approx_equal(post_t1_sim, post_theory, 0.01));
        assert!(approx_equal(post_t2_sim, post_theory, 0.01));
    }

    #[test]
    fn test_small_annotations(){
        let state = vec![false, false, false, false, false];
        let observations = vec![false, true, true, true, false];
        let state_to_observations : Vec<Vec<usize>> = vec![
            vec![0],
            vec![1, 2],
            vec![2, 3],
            vec![3, 4],
            vec![4]];

        let p = 0.5;
        let beta = 0.1;
        let alpha = 0.05;

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta
        );

        let mut state = MgsaState::new(state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Frequency = algorithm.sample(&mut state);

        println!{"{:?}", measure.frequencies}
    }
}
