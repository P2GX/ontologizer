use crate::core::AnnotationIndex;
use ontolius::ontology::OntologyTerms;
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::term::MinimalTerm;
use serde::{Deserialize, Serialize, Serializer};

pub trait Measure {
    /// Returns an iterator over the score for each term.
    fn score(&self) -> f64;

    fn diagnostics(&self) -> Option<String>;
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
    pub fn from_measures<M: Measure>(
        measures: &Vec<M>,
        ontology: &FullCsrOntology,
        annotation_index: &AnnotationIndex,
        observed_genes: &Vec<bool>,
    ) -> Self {
        let term_map = annotation_index.get_terms();
        let gene_map = annotation_index.get_genes();
        let terms_to_genes = annotation_index.get_terms_to_genes(true);

        let mut items = Vec::new();

        for ((measure, term_id), gene_indices) in measures
            .iter()
            .zip(term_map.iter())
            .zip(terms_to_genes.iter())
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
                score: measure.score(),
                associated_genes: gene_symbols,
                diagnostics: measure.diagnostics(),
            })
        }

        EnrichmentResult { items }
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
