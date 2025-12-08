use crate::core::{AnnotationIndex, GeneSymbol};
use ontolius::TermId;
use rand::Rng;
use std::collections::{HashMap, HashSet};

pub fn init_term_states(terms: &HashSet<TermId>, q: f64) -> HashMap<TermId, bool> {
    let mut rng = rand::rng();
    terms
        .iter()
        .map(|t| (t.clone(), rng.random_range(0.0..1.0) < q))
        .collect()
}

pub fn init_hidden_states(genes: &HashSet<GeneSymbol>) -> HashMap<GeneSymbol, bool> {
    genes.iter().map(|g| (g.clone(), false)).collect()
}

pub fn update_hidden_state(
    hidden: &mut HashMap<GeneSymbol, bool>,
    terms: &HashMap<TermId, bool>,
    ann: &AnnotationIndex,
) {
    // reset all genes to false
    for v in hidden.values_mut() {
        *v = false;
    }
    // activate genes annotated to any ON term
    for (term, &on) in terms {
        if on {
            if let Some(genes) = ann.terms_to_genes.get(term) {
                for g in genes {
                    *hidden.get_mut(g).unwrap() = true;
                }
            }
        }
    }
}

pub fn init_obs_states(
    genes: &HashSet<GeneSymbol>,
    obs_genes: &HashSet<GeneSymbol>,
) -> HashMap<GeneSymbol, bool> {
    genes
        .iter()
        .map(|g| (g.clone(), obs_genes.contains(g)))
        .collect()
}
