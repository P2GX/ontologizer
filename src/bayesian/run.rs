#[cfg(test)]
mod test {
    use crate::bayesian::algorithm::{Algorithm, MetropolisHasting};
    use crate::bayesian::model::OrModel;
    use crate::bayesian::proposer::{UniformProposer, UniformToggleProposer};
    use crate::bayesian::recorder::Probability;
    use std::process;

    use crate::bayesian::state::MgsaState;
    use crate::core::{AnnotationIndex, load_gene_set};

    use crate::core::result::EnrichmentResult;
    use csv::Writer;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use ontolius::ontology::csr::FullCsrOntology;

    fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn test_network(
        raw_state: Vec<bool>,
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

        let mut state = MgsaState::new(raw_state);
        let proposer = UniformToggleProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Probability = algorithm.sample(&mut state);
        println! {"{:?}", measure.probabilities}

        for (&sim, theo) in measure.probabilities.iter().zip(posterior) {
            assert!(approx_equal(sim, theo, 0.05));
        }
    }

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

        // MGSA Parameter
        let p = 0.01;
        let alpha = 0.05;
        let beta = 0.10;

        let terms_to_genes = annotation_index.get_terms_to_genes(true);
        let n_genes = annotation_index.genes().len();
        let observed_genes: Vec<bool> = (0..n_genes)
            .map(|i| study_genes.contains(annotation_index.get_index_gene(i)))
            .collect();

        let model = OrModel::new(
            terms_to_genes.clone(),
            observed_genes.clone(),
            p,
            alpha,
            beta,
        );

        let terms = vec![false; terms_to_genes.len()];
        let mut state = MgsaState::new(terms);
        let proposer = UniformToggleProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 50_000_000, 1_000_000);

        let measure: Probability = algorithm.sample(&mut state);

        // Create the Result (Eagerly resolves all strings)
        let mut result = EnrichmentResult::from_measure(
            &measure,
            &ontology,
            annotation_index.terms(),
            annotation_index.genes(),
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
        let terms_to_genes = annotation_index.get_terms_to_genes(true);

        let n_genes = annotation_index.genes().len();
        let observed_genes: Vec<bool> = (0..n_genes)
            .map(|i| study_genes.contains(annotation_index.get_index_gene(i)))
            .collect();

        let p = 0.01;
        let alpha = 0.05;
        let beta = 0.10;
        let model = OrModel::new(
            terms_to_genes.clone(),
            observed_genes.clone(),
            p,
            alpha,
            beta,
        );

        let terms = vec![false; terms_to_genes.len()];
        let mut state = MgsaState::new(terms);
        let proposer = UniformToggleProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 5_000_000, 1_000_000);

        let measure: Probability = algorithm.sample(&mut state);

        let mut result = EnrichmentResult::from_measure(
            &measure,
            &ontology,
            annotation_index.terms(),
            annotation_index.genes(),
            &observed_genes,
            &terms_to_genes,
            p,
            alpha,
            beta,
        );

        // Optional: Sort
        result.sort_by_score(true); // descending for probability

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GO0090717/results.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }

    #[test]
    fn test_background_noise() {
        fn log_sum_exp(vals: &[f64]) -> f64 {
            let max_val = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let sum = vals.iter().map(|v| (v - max_val).exp()).sum::<f64>();
            max_val + sum.ln()
        }

        let n = 500;

        let p = 0.01;
        let alpha = 0.05;
        let beta = 0.10;

        let state = vec![false; n];
        let mut state_to_observations: Vec<Vec<usize>> = (0..n).map(|_| (0..n).collect()).collect();

        // Observations
        let mut observations = vec![false; n];

        observations[0] = true;

        state_to_observations[0] = vec![0];

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 500_000, 100_000);

        let measure: Probability = algorithm.sample(&mut state);

        let n_background = n - 1; // 49

        // Log-probabilities helper
        let ln = |x: f64| x.ln();

        // Priors (Log space)
        // P(t1=1) and P(t1=0)
        let lp_t1_1 = ln(p);
        let lp_t1_0 = ln(1.0 - p);

        // P(S=0) = (1-p)^49
        let lp_s_0 = (n_background as f64) * ln(1.0 - p);
        // P(S=1) = 1 - (1-p)^49. computed in linear space for precision then logged
        let p_s_1 = 1.0 - (1.0 - p).powi(n_background as i32);
        let lp_s_1 = ln(p_s_1);

        // Likelihoods (Log space)
        // Common factors
        let ln_alpha = ln(alpha);
        let ln_not_alpha = ln(1.0 - alpha);
        let ln_beta = ln(beta);
        let ln_not_beta = ln(1.0 - beta); // 1-beta is sensitivity (0.9)

        // L_00: t1=0, S=0. h1=0, h_rest=0.
        // o1=1 (FP), o_rest=0 (TN)
        let ll_00 = ln_alpha + (n_background as f64) * ln_not_alpha;

        // L_10: t1=1, S=0. h1=1, h_rest=0.
        // o1=1 (TP), o_rest=0 (TN)
        let ll_10 = ln_not_beta + (n_background as f64) * ln_not_alpha;

        // L_01: t1=0, S=1. h1=1, h_rest=1.
        // o1=1 (TP), o_rest=0 (FN)
        let ll_01 = ln_not_beta + (n_background as f64) * ln_beta;

        // L_11: t1=1, S=1. h1=1, h_rest=1.
        // Same likelihood as L_01
        let ll_11 = ll_01;

        // Joint Log Probabilities: Likelihood + Prior
        let lj_00 = ll_00 + lp_t1_0 + lp_s_0;
        let lj_10 = ll_10 + lp_t1_1 + lp_s_0;
        let lj_01 = ll_01 + lp_t1_0 + lp_s_1;
        let lj_11 = ll_11 + lp_t1_1 + lp_s_1;

        // Log-Sum-Exp to get Evidence E
        let log_evidence = log_sum_exp(&[lj_00, lj_10, lj_01, lj_11]);

        // Posteriors
        // P(t1=1 | o) = (J_10 + J_11) / E
        let log_prob_t1 = log_sum_exp(&[lj_10, lj_11]) - log_evidence;
        let prob_t1 = log_prob_t1.exp();

        // P(S=1 | o)
        let log_prob_s = log_sum_exp(&[lj_01, lj_11]) - log_evidence;
        let prob_s = log_prob_s.exp();

        // P(tk=1 | o) = P(S=1 | o) * p / P(S=1)
        let prob_tk = prob_s * (p / p_s_1);

        println!("{:?} {:?}", prob_t1, measure.probabilities[0]);
        println!("{:?} {:?}", prob_tk, measure.probabilities[1]);
        assert!(approx_equal(prob_t1, measure.probabilities[0], 0.05));
        assert!(approx_equal(prob_tk, measure.probabilities[1], 0.05));
    }

    #[test]
    fn test_term_broadness() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let raw_state = vec![false, false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![
            vec![0, 1, 2],
            vec![0, 1, 2, 3, 4],
            vec![0, 1, 2, 3, 4, 5, 6],
        ];

        // Observations
        let observations = vec![true, true, true, false, false, false, false];

        let model = OrModel::new(
            state_to_observations.clone(),
            observations.clone(),
            p,
            alpha,
            beta,
        );

        let mut state = MgsaState::new(raw_state);
        let proposer = UniformProposer::new();
        let mut algorithm = MetropolisHasting::new(model, proposer, 10_000, 0);

        let measure: Probability = algorithm.sample(&mut state);
        println! {"{:?}", measure.probabilities}
        assert!(measure.probabilities[0] > measure.probabilities[1]);
        assert!(measure.probabilities[1] > measure.probabilities[2]);
    }

    #[test]
    fn test_one_term_one_gene() {
        let p = 0.05;
        let alpha = 0.05;
        let beta = 0.1;

        let raw_state = vec![false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0]];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior = vec![p * (1. - beta) / (alpha * (1. - p) + (1. - beta) * p)];

        test_network(
            raw_state.clone(),
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
            raw_state.clone(),
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
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0], vec![0]];

        // Active Observation O=(1)
        let observations = vec![true];
        let posterior_t1 =
            p * (1. - beta) / (alpha * (1. - p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));
        let posterior_t2 =
            p * (1. - beta) / (alpha * (1. - p).powf(2.) + (1. - beta) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false];
        let posterior_t1 =
            p * beta / ((1. - alpha) * (1. - p).powf(2.) + beta * (1. - (1. - p).powf(2.)));
        let posterior_t2 =
            p * beta / ((1. - alpha) * (1. - p).powf(2.) + beta * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
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
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
        let state_to_observations: Vec<Vec<usize>> = vec![vec![0, 1], vec![0, 1]];

        // Active Observation O=(1)
        let observations = vec![true, true];
        let posterior_t1 = p * (1. - beta).powf(2.)
            / (alpha.powf(2.) * (1. - p).powf(2.)
                + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * (1. - beta).powf(2.)
            / (alpha.powf(2.) * (1. - p).powf(2.)
                + (1. - beta).powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            posterior,
        );

        // Inactive Observation O=(0)
        let observations = vec![false, false];
        let posterior_t1 = p * beta.powf(2.)
            / ((1. - alpha).powf(2.) * (1. - p).powf(2.)
                + beta.powf(2.) * (1. - (1. - p).powf(2.)));
        let posterior_t2 = p * beta.powf(2.)
            / ((1. - alpha).powf(2.) * (1. - p).powf(2.)
                + beta.powf(2.) * (1. - (1. - p).powf(2.)));

        let posterior = vec![posterior_t1, posterior_t2];

        test_network(
            raw_state.clone(),
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
        let p: f64 = 0.05;
        let alpha: f64 = 0.05;
        let beta: f64 = 0.1;

        let raw_state = vec![false, false];
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
            raw_state.clone(),
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
            raw_state.clone(),
            &observations,
            &state_to_observations,
            p,
            alpha,
            beta,
            vec![post_t0, post_t1],
        );
    }
}
