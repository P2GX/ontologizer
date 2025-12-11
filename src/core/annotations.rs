use crate::core::geneset::{GeneSet, GeneSymbol};
use oboannotation::{
    go::{GoAnnotations, GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use ontolius::{
    TermId,
    ontology::{HierarchyWalks, csr::FullCsrOntology},
};
use std::collections::{HashMap, HashSet};

// Contains all GO annotations and provides methods to build annotation maps specific to the study and population sets.
pub struct AnnotationIndex {
    pub annotations: GoAnnotations, // All GO annotations loaded from the GAF file
    pub terms_to_genes: HashMap<TermId, HashSet<GeneSymbol>>, // all GO terms mapped to their annotated genes
    pub genes_to_terms: HashMap<GeneSymbol, HashSet<TermId>>, // all genes mapped to their annotated GO terms
}

impl AnnotationIndex {
    pub fn new(gaf_path: &str, go: &FullCsrOntology) -> Self {
        // Load annotations from GAF file
        let loader = GoGafAnnotationLoader;
        let annotations = loader
            .load_from_path(gaf_path)
            .expect("Could not load GAF file");

        let genes_to_terms: HashMap<GeneSymbol, HashSet<TermId>> = get_annotation_map(&annotations)
            .into_iter()
            .filter_map(|(k, v)| match GeneSymbol::new(&k) {
                Some(sym) => Some((sym, v)),
                None => {
                    eprintln!("Could not parse gene symbol: {k}");
                    None
                }
            })
            .collect();

        let mut terms_to_genes: HashMap<TermId, HashSet<GeneSymbol>> = HashMap::new();

        for (gene, terms) in &genes_to_terms {
            for term in terms {
                for ancestor in go.iter_term_and_ancestor_ids(term) {
                    terms_to_genes
                        .entry(ancestor.clone())
                        .or_default()
                        .insert(gene.clone());
                }
            }
        }

        Self {
            annotations,
            genes_to_terms,
            terms_to_genes,
        }
    }

    pub fn terms(&self) -> HashSet<TermId> {
        self.terms_to_genes.keys().cloned().collect()
    }

    pub fn terms_for_subset(&self, genes: &GeneSet, ontology: &FullCsrOntology) -> HashSet<TermId> {
        let mut terms_subset = HashSet::new();

        for gene in genes.recognized_genes() {
            for term in self.genes_to_terms.get(gene).unwrap() {
                for ancestor in ontology.iter_term_and_ancestor_ids(term) {
                    terms_subset.insert(ancestor.clone());
                }
            }
        }
        terms_subset
    }

    pub fn terms_to_genes_for_subset(
        &self,
        genes: &GeneSet,
        ontology: &FullCsrOntology,
    ) -> HashMap<TermId, HashSet<GeneSymbol>> {
        let mut terms_to_genes: HashMap<TermId, HashSet<GeneSymbol>> = HashMap::new();

        for gene in genes.recognized_genes() {
            let mut terms = HashSet::new();
            for term in self.genes_to_terms.get(gene).unwrap() {
                terms.extend(ontology.iter_term_and_ancestor_ids(term).cloned()); // true-path rule
            }
            for term in terms {
                terms_to_genes
                    .entry(term.clone())
                    .or_default()
                    .insert(gene.clone());
            }
        }
        terms_to_genes
    }

    pub fn term_counts_for_subset(
        &self,
        genes: &GeneSet,
        ontology: &FullCsrOntology,
    ) -> HashMap<TermId, usize> {
        let mut counts = HashMap::new();
        for gene in genes.recognized_genes() {
            let mut terms = HashSet::new();
            for term in self.genes_to_terms.get(gene).unwrap() {
                terms.extend(ontology.iter_term_and_ancestor_ids(term).cloned()); // true-path rule
            }
            for term in terms {
                *counts.entry(term).or_insert(0) += 1;
            }
        }
        counts
    }

    pub fn get_genes_to_terms(&self) -> &HashMap<GeneSymbol, HashSet<TermId>> {
        &self.genes_to_terms
    }

    pub fn get_terms_to_genes(&self) -> &HashMap<TermId, HashSet<GeneSymbol>> {
        &self.terms_to_genes
    }

    pub fn get_annotations(&self) -> &GoAnnotations {
        &self.annotations
    }
}
