use indexmap::IndexSet;
use ontolius::TermId;
use ontolius::ontology::OntologyTerms;
use ontolius::ontology::csr::FullCsrOntology;
use serde::{Deserialize, Serialize, Serializer};
use std::collections::HashMap;
use ontolius::term::MinimalTerm;

pub trait Measure {
    /// Returns an iterator over the score for each term.
    fn scores(&self) -> impl Iterator<Item = f64>;

    fn diagnostics(&self) -> impl Iterator<Item = Option<String>>;

    fn get_score(&self, i: usize) -> f64;

    fn get_diagnostics(&self, i: usize) -> Option<f64>;
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

    /// Arbitrary additional info
    #[serde(rename = "Diagnostics", skip_serializing_if = "Option::is_none")]
    pub diagnostics: Option<String>,
}

fn serialize_genes<S>(genes: &Vec<String>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let joined = genes.join(", ");
    serializer.serialize_str(&joined)
}

#[derive(Serialize)]
pub struct EnrichmentResult {
    pub items: Vec<EnrichmentItem>,
}

impl EnrichmentResult {
    pub fn from_measure<M: Measure>(
        measures: &M,
        ontology: &FullCsrOntology,
        term_map: &IndexSet<TermId>,
        gene_map: &IndexSet<String>,
        observed_genes: &Vec<bool>,
        terms_to_genes: &Vec<Vec<usize>>,
    ) -> Self {
        let mut items = Vec::new();

        for (((score, term_id), gene_indices), diagnostics) in measures
            .scores()
            .zip(term_map.iter())
            .zip(terms_to_genes.iter())
            .zip(measures.diagnostics())
        {
            let label = ontology
                .term_by_id(term_id)
                .map(|t| t.name().to_string())
                .unwrap_or_else(|| "Unknown Term".to_string());

            let gene_symbols: Vec<String> = gene_indices
                .iter()
                .filter(|&idx| observed_genes[*idx])
                .map(|&idx| gene_map[idx].to_string())
                .collect();

            items.push(EnrichmentItem {
                id: term_id.to_string(),
                label,
                score,
                associated_genes: gene_symbols,
                diagnostics,
            })
        }
        
        EnrichmentResult {
            items,
        }
    }


    pub fn sort_by_score(&mut self, descending: bool) {
        if descending {
            self.items
                .sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        } else {
            self.items
                .sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap())
        }
    }
}