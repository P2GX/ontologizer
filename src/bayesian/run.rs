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

    fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn test_network(
        observations: &Vec<bool>,
        state_to_observations: &Vec<Vec<usize>>,
        p: f64,
        alpha: f64,
        beta: f64,
        posterior: Vec<f64>,
    ) {
        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let terms = model.heuristic_start();
        let mut state = MgsaState::new(terms);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Frequency = algorithm.sample(&mut state);
        println! {"{:?}", measure.frequencies}

        for (&sim, theo) in measure.frequencies.iter().zip(posterior) {
            assert!(approx_equal(sim, theo, 0.01));
        }
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
    fn test_one_term_one_gene() {
        let p = 0.5;
        let alpha = 0.05;
        let beta = 0.1;

        let state = vec![false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0]];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior = vec![p * (1. - beta) / (alpha * (1. - p) + (1. - beta) * p)];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false];
        let posterior = vec![p * beta / ((1. - alpha) * (1. - p) + beta * p)];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_two_term_one_gene() {
        let p : f64 = 0.5;
        let alpha : f64 = 0.05;
        let beta  : f64 = 0.1;

        let state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0],
            vec![0]
        ];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior_t1 = p * (1. - beta) / (alpha*(1.- p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * (1. - beta) / (alpha*(1.- p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false];
        let posterior_t1 = p * beta / ((1.-alpha)*(1.- p).powf(2.) + beta * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * beta / ((1.-alpha)*(1.- p).powf(2.) + beta * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_full_two_term_two_gene() {
        let p : f64 = 0.5;
        let alpha : f64 = 0.05;
        let beta  : f64 = 0.1;

        let state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0, 1],
            vec![0, 1]
        ];

        // Active Observation O=(1)
        let observations = vec![true, true];
        let posterior_t1 = p * (1. - beta).powf(2.) / (alpha.powf(2.)*(1.- p).powf(2.) + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * (1. - beta).powf(2.) / (alpha.powf(2.)*(1.- p).powf(2.) + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false, false];
        let posterior_t1 = p * beta.powf(2.) / ((1.-alpha).powf(2.)*(1.- p).powf(2.) + beta.powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * beta.powf(2.) / ((1.-alpha).powf(2.)*(1.- p).powf(2.) + beta.powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );
    }

    #[test]
    fn test_asymmetric_two_term_two_gene() {
        let p: f64 = 0.5;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0],    // Term 0
            vec![0, 1], // Term 1
        ];

        let observations = vec![true, false];

        let lh_00 = alpha * (1.0 - alpha);
        let lh_10 = (1.0 - beta) * (1.0 - alpha);
        let lh_01 = (1.0 - beta) * beta;
        let lh_11 = (1.0 - beta) * beta;

        let joint_00 = lh_00 * (1.0 - p).powf(2.0);
        let joint_10 = lh_10 * p * (1.0 - p);
        let joint_01 = lh_01 * p * (1.0 - p);
        let joint_11 = lh_11 * p.powf(2.0);

        let z = joint_00 + joint_10 + joint_01 + joint_11;

        // Marginalize for Posteriors
        let post_t0 = (joint_10 + joint_11) / z;
        let post_t1 = (joint_01 + joint_11) / z;

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );


        let observations = vec![false, true];

        let lh_00 = (1.0 - alpha) * alpha;
        let lh_10 = beta * alpha;
        let lh_01 = beta * (1.0 - beta);
        let lh_11 = beta * (1.0 - beta);

        let joint_00 = lh_00 * (1.0 - p).powf(2.0);
        let joint_10 = lh_10 * p * (1.0 - p);
        let joint_01 = lh_01 * p * (1.0 - p);
        let joint_11 = lh_11 * p.powf(2.0);

        let z = joint_00 + joint_10 + joint_01 + joint_11;

        let post_t0 = (joint_10 + joint_11) / z;
        let post_t1 = (joint_01 + joint_11) / z;

        test_network(
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );
    }
}
