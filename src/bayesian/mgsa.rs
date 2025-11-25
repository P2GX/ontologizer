use crate::bayesian::mcmc::Sampler;
use crate::bayesian::probs::*;
use crate::bayesian::types::{MgsaConfig, MgsaParams, MgsaResult, Problem, State};
use crate::bayesian::utils::{
    init_hidden_states, init_obs_states, init_term_states, update_hidden_state,
};
use rand::Rng;

pub fn run_mgsa(problem: &Problem, params: &MgsaParams, cfg: &MgsaConfig) -> MgsaResult {
    // todo!(should allow to restrict to a subset of terms)
    let obs_genes = problem.observed_genes();
    let all_genes = problem.all_genes();

    // todo!(should allow to restrict to a subset of terms - determined by genes)
    let all_terms = problem.all_terms();

    let mut result = MgsaResult::new(&all_terms);

    // initial state: all terms off, hidden inferred from that
    let mut init_term_states = init_term_states(&all_terms, params.q);
    let mut init_hidden_states = init_hidden_states(&all_genes);
    let obs_states = init_obs_states(&all_genes, &obs_genes);

    update_hidden_state(&mut init_hidden_states, &init_term_states, &problem.ann);

    let mut state = State {
        terms: init_term_states,
        hidden: init_hidden_states,
    };

    let mut state_p = state.clone();

    let mut sampler = Sampler::new(); // todo!(allow to configure with parameters)
    let false_pos = params.alpha; // false positive rate (O_i = 1 | H_i = 0)
    let false_neg = params.beta; // false negative rate (O_i = 0 | H_i = 1)
    let pri_prob = params.q; // a priori probability p(T_i = 1)

    // precompute current score
    // todo!(in the entire code below: function signature should be optimized to utilize types.rs)
    let mut log_post_t = log_posterior_unnorm(
        false_pos,
        false_neg,
        pri_prob,
        &obs_states,
        &state.hidden,
        &state.terms,
    );
    let mut neigh_t = neighborhood_size(&state.terms) as f64;

    let mut rng = rand::rng();
    for step in 0..cfg.steps {
        if step % 1_000 == 0 {
            println!("MGSA sampling step {}/{}", step, cfg.steps);
        }

        // propose terms to be altered
        let term_move = sampler.draw_move(&state.terms);
        // rebuild hidden from proposed terms
        sampler.apply_move(&mut state_p.terms, &term_move);

        // build Metropolis–Hastings acceptance probability for the proposed state
        let log_post_p = log_posterior_unnorm(
            false_pos,
            false_neg,
            pri_prob,
            &obs_states,
            &state_p.hidden,
            &state_p.terms,
        );
        // proposal probability q(t^p | t^t) = 1 / N^t
        let neigh_p = neighborhood_size(&state_p.terms) as f64;

        // acceptance ratio in log-space:
        // log a = (log_post_prop - log_post_curr) + log(neigh_curr/neigh_prop)
        // prob = proposed state, curr = current state
        let log_a = (log_post_t - log_post_p) + neigh_t.ln() - neigh_p.ln();

        let x: f64 = rng.random_range(0.0..1.0);

        if x.ln() < log_a {
            // accept the proposed state
            log_post_t = log_post_p;
            neigh_t = neigh_p;
        } else {
            // reject the proposed state
            sampler.revert_move(&mut state_p.terms, &term_move);
        }

        if step >= cfg.burn_in {
            result.update(&state);
        }
    }

    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::bayesian::types::{MgsaParams, Problem};
    use crate::core::AnnotationIndex;
    use crate::core::Ontologizer;
    use crate::core::{load_gene_set, separate_gene_set};

    #[test]
    fn test_mgsa() {
        println!("Testing MGSA");

        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let study_set_path = "tests/data/study.txt";
        let pop_set_path = "tests/data/population.txt";

        let ontologizer = Ontologizer::new(go_path);
        let onto = ontologizer.ontology();
        let mut ann = AnnotationIndex::new(gaf_path, onto);

        // Load the population and study gene sets
        let pop_gene_symbols =
            load_gene_set(pop_set_path).expect("Failed to parse population gene set");
        let pop_gene_set = separate_gene_set(&ann.get_annotations(), pop_gene_symbols);
        let uni_genes = pop_gene_set.recognized_genes().clone();

        let study_gene_symbols =
            load_gene_set(study_set_path).expect("Failed to parse study gene set");
        let study_gene_set = separate_gene_set(&ann.get_annotations(), study_gene_symbols);
        let obs_genes = study_gene_set.recognized_genes().clone();

        // Fix term universe and order
        // todo!(if a population is provided, limit the terms to those annotated to the population)
        let terms = ann.terms();

        // Bind the problem view
        let problem = Problem {
            ann: &ann,
            genes: obs_genes,
        };

        // Configure fixed parameters and MCMC
        let params = MgsaParams {
            alpha: 0.10, // false positive rate
            beta: 0.05,  // false negative rate
            q: 0.1,      // prior prob term is on
        };

        let cfg = MgsaConfig {
            steps: 20_000,
            burn_in: 10_000,
        };

        // Run sampler
        println!("Running MGSA sampler");
        let out = run_mgsa(&problem, &params, &cfg);

        out.write_tsv("tests/data/enrichment_results_mgsa.tsv")
            .expect("Failed to write results");
    }
}
