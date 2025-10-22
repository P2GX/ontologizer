use std::collections::{HashMap, HashSet};

use crate::core::geneset::{GeneSet, GeneSymbol};
use oboannotation::{
    go::{GoAnnotations, GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use ontolius::{
    TermId,
    ontology::{HierarchyWalks, csr::FullCsrOntology},
};

// Contains all GO annotations and provides methods to build annotation maps specific to the study and population sets.
pub struct AnnotationIndex {
    pub annotations: GoAnnotations, // All GO annotations loaded from the GAF file
    pub genes_to_terms: HashMap<GeneSymbol, HashSet<TermId>>,
    pub terms_to_genes: HashMap<TermId, HashSet<GeneSymbol>>,
    pub population_term_counts: HashMap<TermId, usize>, // Counts of all terms annotated in the population set
    pub study_term_counts: HashMap<TermId, usize>, // Counts of all terms annotated in the study set
}

impl AnnotationIndex {
    pub fn new(gaf_path: &str) -> Self {
        // Load annotations from GAF file
        let loader = GoGafAnnotationLoader;
        let annotations = loader
            .load_from_path(gaf_path)
            .expect("Could not load GAF file");


        let genes_to_terms: HashMap<GeneSymbol, HashSet<TermId>> = get_annotation_map(&annotations)
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

        let terms_to_genes : HashMap<TermId, HashSet<GeneSymbol>> = invert_annotation_map(&genes_to_terms);
        
        Self {
            annotations,
            genes_to_terms,
            terms_to_genes,
            population_term_counts: HashMap::new(),
            study_term_counts: HashMap::new(),
        }
    }

    // Builds a hashmap of GO terms annotated to a gene set and their counts.
    pub fn compute_term_counts(&mut self, gene_set: &GeneSet, go: &FullCsrOntology) {
        let mut counts = HashMap::new();
        for gene in gene_set.recognized_genes() {
            let mut terms = HashSet::new();
            for term in self.genes_to_terms.get(gene).unwrap() {
                terms.extend(go.iter_term_and_ancestor_ids(term).cloned()); // true-path rule
            }
            for term in terms {
                *counts.entry(term).or_insert(0) += 1;
            }
        }
        self.study_term_counts = counts;
    }


    // Builds a map of all GO terms that are annotated to genes in the population set and their associated genes.
    // We iterate through all genes in the population set and check which terms they are directly annotated to.
    // For each of these directly annotated terms and all of their ancestors, we add the term and the gene to the map.
    // This means for each term the length of the gene set is mt
    pub fn build_terms_to_genes(&mut self, pop_set: &GeneSet, go: &FullCsrOntology) {
        let mut terms_to_genes: HashMap<TermId, HashSet<GeneSymbol>> = HashMap::new();

        for gene in pop_set.recognized_genes() {
            if let Some(direct_terms) = self.genes_to_terms.get(gene) {
                let mut seen_terms = HashSet::new();
                for term in direct_terms {
                    for ancestor in go.iter_term_and_ancestor_ids(term) {
                        // Verhindert mehrfaches Einfügen für denselben Term
                        if seen_terms.insert(ancestor.clone()) {
                            terms_to_genes
                                .entry(ancestor.clone())
                                .or_insert_with(HashSet::new)
                                .insert(gene.clone());
                        }
                    }
                }
            }
        }
        self.terms_to_genes = terms_to_genes;
    }

    pub fn study_annotations(&self) -> &HashMap<TermId, usize> {
        &self.study_term_counts
    }

    pub fn annotation_map(&self) -> &HashMap<GeneSymbol, HashSet<TermId>> {
        &self.genes_to_terms
    }

    pub fn term_genes_map(&self) -> &HashMap<TermId, HashSet<GeneSymbol>> {
        &self.terms_to_genes
    }

    pub fn annotations(&self) -> &GoAnnotations {
        &self.annotations
    }
}

 fn invert_annotation_map(map: &HashMap<GeneSymbol, HashSet<TermId>>) -> HashMap<TermId, HashSet<GeneSymbol>> {
     let mut inv_map : HashMap<TermId, HashSet<GeneSymbol>> = HashMap::new();
     for (gene, terms) in map {
         for term in terms {
             inv_map.entry(term.clone()).or_default().insert(gene.clone());
         }
     }
     inv_map
}

/* pub fn get_inverted_annotation_map(annotations: &GoAnnotations) -> HashMap<String, HashSet<TermId>> {
    let mut inv_annot_map = HashMap::new();
    for ann in &annotations.annotations {
        if ann.is_negated() {
            continue;
        }
        let symbol = ann.gene_product_symbol.clone();
        let tid = ann.gene_ontology_id.clone();
        // Insert into HashMap, creating a new HashSet if necessary
        inv_annot_map
            .entry(tid)
            .or_insert_with(|| HashSet::new())
            .insert(symbol); //
    }
    inv_annot_map
} */