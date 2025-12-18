use crate::core::GeneSymbol;
use indexmap::IndexSet;
use ontolius::TermId;
use ontolius::ontology::OntologyTerms;
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::term::MinimalTerm;
use serde::{Deserialize, Serialize, Serializer};
use std::collections::HashMap;

pub trait AnalysisResult {
    /// Returns the main table of results.
    fn items(&self) -> &[EnrichmentItem];

    /// Returns metadata about the run (e.g., "p=0.5", "alpha=0.1", "Correction=FDR").
    fn parameters(&self) -> HashMap<String, String>;

    /// Helper to sort results by score (descending for Prob, ascending for P-val).
    fn sort_by_score(&mut self, descending: bool);
}

pub trait Measure {
    /// Returns an iterator over the score for each term.
    fn scores(&self) -> impl Iterator<Item = f64>;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnrichmentItem {
    #[serde(rename = "ID")]
    pub id: String,

    #[serde(rename = "Label")]
    pub label: String,

    /// The primary metric: Posterior Probability (Bayesian) or P-Value (Frequentist)
    #[serde(rename = "Score")]
    pub score: f64,

    /// Genes from the study set annotated to this term
    #[serde(rename = "Associated Genes", serialize_with = "serialize_genes")]
    pub associated_genes: Vec<String>,
}

fn serialize_genes<S>(genes: &Vec<String>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let joined = genes.join(", ");
    serializer.serialize_str(&joined)
}

#[derive(Serialize)]
pub struct BayesianResult {
    parameters: HashMap<String, String>,
    items: Vec<EnrichmentItem>,
}

impl BayesianResult {
    pub fn from_counts<M: Measure>(
        measures: &M,
        ontology: &FullCsrOntology,
        term_map: &IndexSet<TermId>,
        gene_map: &IndexSet<GeneSymbol>,
        terms_to_genes: &Vec<Vec<usize>>,
        // Metadata
        p: f64,
        alpha: f64,
        beta: f64,
    ) -> Self {
        let mut items = Vec::new();

        for ((measure, term_id), gene_indices) in measures
            .scores()
            .zip(term_map.iter())
            .zip(terms_to_genes.iter())
        {
            let label = ontology
                .term_by_id(term_id)
                .map(|t| t.name().to_string())
                .unwrap_or_else(|| "Unknown Term".to_string());

            let gene_symbols: Vec<String> = gene_indices
                .iter()
                .map(|&idx| gene_map[idx].to_string())
                .collect();

            items.push(EnrichmentItem {
                id: term_id.to_string(),
                label,
                score: measure,
                associated_genes: gene_symbols,
            })
        }

        let mut params = HashMap::new();
        params.insert("p".to_string(), p.to_string());
        params.insert("alpha".to_string(), alpha.to_string());
        params.insert("beta".to_string(), beta.to_string());

        BayesianResult {
            parameters: params,
            items,
        }
    }
}

impl AnalysisResult for BayesianResult {
    fn items(&self) -> &[EnrichmentItem] {
        &self.items
    }

    fn parameters(&self) -> HashMap<String, String> {
        self.parameters.clone()
    }

    fn sort_by_score(&mut self, descending: bool) {
        if descending {
            self.items
                .sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        } else {
            self.items
                .sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
        }
    }
}
