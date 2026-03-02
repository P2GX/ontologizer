use crate::core::AnnotationIndex;
use csv::Writer;
use ontolius::TermId;
use ontolius::common::go::{BIOLOGICAL_PROCESS, CELLULAR_COMPONENT, MOLECULAR_FUNCTION};
use ontolius::ontology::csr::FullCsrOntology;
use ontolius::ontology::{HierarchyQueries, OntologyTerms};
use ontolius::term::MinimalTerm;
use serde::{Deserialize, Serialize, Serializer};
use std::path::Path;
use std::slice::Iter;

pub trait Measure {
    /// Returns an iterator over the score for each term.
    fn score(&self) -> f64;

    fn diagnostics(&self) -> Option<String>;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnrichmentItem {
    #[serde(rename = "Id")]
    pub id: String,

    #[serde(rename = "Label")]
    pub label: String,

    #[serde(rename = "Aspect")]
    pub aspect: String,

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

            let aspect = get_term_aspect(ontology, term_id);

            let gene_symbols: Vec<String> = gene_indices
                .iter()
                .filter(|&idx| observed_genes[*idx])
                .map(|&idx| gene_map[idx].to_string())
                .collect();

            items.push(EnrichmentItem {
                id: term_id.to_string(),
                label,
                aspect,
                score: measure.score(),
                associated_genes: gene_symbols,
                diagnostics: measure.diagnostics(),
            })
        }

        EnrichmentResult { items }
    }

    pub fn iter(&self) -> Iter<'_, EnrichmentItem> {
        self.items.iter()
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

    /// Saves the result to a CSV file.
    ///
    /// This writes the global parameters as commented lines (header)
    /// before writing the standard CSV table.
    pub fn save_to_csv<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let mut wtr = Writer::from_path("parameters.csv")?;
        for item in &self.items {
            wtr.serialize(item)?;
        }

        let mut wtr = Writer::from_path("terms.csv")?;
        for item in &self.items {
            wtr.serialize(item)?;
        }
        Ok(())
    }
}

fn get_term_aspect(go: &FullCsrOntology, term: &TermId) -> String {
    if go.is_equal_or_descendant_of(term, &BIOLOGICAL_PROCESS) {
        "BP".to_string()
    } else if go.is_equal_or_descendant_of(term, &MOLECULAR_FUNCTION) {
        "MF".to_string()
    } else if go.is_equal_or_descendant_of(term, &CELLULAR_COMPONENT) {
        "CC".to_string()
    } else {
        "Unknown".to_string()
    }
}
