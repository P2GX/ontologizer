use std::collections::{HashMap, HashSet};

use crate::geneset::{GeneSet, GeneSymbol};
use oboannotation::{
    go::{GoAnnotations, GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use ontolius::{
    TermId,
    ontology::{HierarchyWalks, csr::FullCsrOntology},
};
use std::time::Instant;

// Contains all GO annotations and provides methods to build annotation maps specific to the study and population sets.
pub struct AnnotationContainer {
    pub annotations: GoAnnotations, // Contains all GO annotations loaded from the GAF file
    pub annotation_map: HashMap<GeneSymbol, HashSet<TermId>>, // Maps genes to directly annotated GO terms
    study_annotations: HashMap<TermId, usize>, // Maps GO terms annotated in the study set to their counts in the study set (nt)
    term_genes_map: HashMap<TermId, HashSet<GeneSymbol>>, // Maps GO terms annotated in the population set to their associated genes in the population set
}

impl AnnotationContainer {
    pub fn new(gaf_path: &str) -> Self {
        let loader = GoGafAnnotationLoader;
        let annotations = loader
            .load_from_path(gaf_path)
            .expect("Could not load GAF file");

        let annotation_map: HashMap<GeneSymbol, HashSet<TermId>> = get_annotation_map(&annotations)
            .into_iter()
            .filter_map(|(k, v)| {
                match GeneSymbol::new(&k) {
                    Some(sym) => Some((sym, v)),
                    None => {
                        eprintln!("Could not parse gene symbol: {k}");
                        None
                    }
                }
            })
            .collect();

        Self {
            annotations,
            annotation_map: annotation_map,
            study_annotations: HashMap::new(),
            term_genes_map: HashMap::new(),
        }
    }

    // Builds a hashmap of all GO terms annotated in the study set and their counts.
    // This way, we know which terms we actually need to analyze.
    pub fn build_study_annotations(&mut self, study_set: &GeneSet, go: &FullCsrOntology) {
        let mut study_annotated_terms = HashMap::new();
        for gene in study_set.recognized_gene_symbols() {
            let mut unique_terms = HashSet::new();

            for term in self.annotation_map.get(gene).unwrap() {
                unique_terms.extend(go.iter_term_and_ancestor_ids(term).cloned());
            }

            for term in unique_terms {
                *study_annotated_terms.entry(term).or_insert(0) += 1;
            }
        }

        self.study_annotations = study_annotated_terms;
    }

    // Builds a map of all GO terms that are annotated to genes in the population set and their associated genes.
    // We iterate through all genes in the population set and check which terms they are directly annotated to.
    // For each of these directly annotated terms and all of their ancestors, we add the term and the gene to the map.
    // This means for each term the length of the gene set is mt
    pub fn build_term_genes_map(&mut self, pop_set: &GeneSet, go: &FullCsrOntology) {
        let mut term_to_genes: HashMap<TermId, HashSet<GeneSymbol>> = HashMap::new();
        let start = Instant::now();

        for gene in pop_set.recognized_gene_symbols() {
            if let Some(direct_terms) = self.annotation_map.get(gene) {
                let mut seen_terms = HashSet::new();
                for term in direct_terms {
                    for ancestor in go.iter_term_and_ancestor_ids(term) {
                        // Verhindert mehrfaches Einfügen für denselben Term
                        if seen_terms.insert(ancestor.clone()) {
                            term_to_genes
                                .entry(ancestor.clone())
                                .or_insert_with(HashSet::new)
                                .insert(gene.clone());
                        }
                    }
                }
            }
        }
        eprintln!("Building term to genes map took: {:?}", start.elapsed());
        self.term_genes_map = term_to_genes;
    }

    pub fn study_annotations(&self) -> &HashMap<TermId, usize> {
        &self.study_annotations
    }

    pub fn annotation_map(&self) -> &HashMap<GeneSymbol, HashSet<TermId>> {
        &self.annotation_map
    }

    pub fn term_genes_map(&self) -> &HashMap<TermId, HashSet<GeneSymbol>> {
        &self.term_genes_map
    }

    pub fn annotations(&self) -> &GoAnnotations {
        &self.annotations
    }
}
