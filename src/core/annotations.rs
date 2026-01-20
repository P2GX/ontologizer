use oboannotation::go::{GoAnnotations, stats::get_annotation_map};
use ontolius::{
    TermId,
    ontology::{HierarchyWalks, csr::FullCsrOntology},
};

use indexmap::IndexSet;
use std::collections::{HashMap, HashSet};

// Contains all GO annotations and provides methods to build annotation maps specific to the population sets.
pub struct AnnotationIndex {
    pub annotations: GoAnnotations, // All GO annotations loaded from the GAF file

    // Bidirectional mapping
    // Index -> String: self.terms.get_index(i)
    // String -> Index: self.terms.get_index_of(term_id)
    term_map: IndexSet<TermId>,
    gene_map: IndexSet<String>,

    // The Graph represented by integer adjacency matrices
    // gene_to_term[j] = list of term indices gene j is annotated to
    // term_to_gene[i] = list of gene indices annotated to term i

    // sparse / excluding ancestors / implicit
    gene_to_term_sparse: Vec<IndexSet<usize>>,
    term_to_gene_sparse: Vec<IndexSet<usize>>,

    // dense / including ancestors / explicit
    gene_to_term_dense: Vec<IndexSet<usize>>,
    term_to_gene_dense: Vec<IndexSet<usize>>,
}

impl AnnotationIndex {
    /// Creates a new instance of the data structure by processing gene-to-term annotations,
    /// applying optional population gene filtering, and constructing all relevant mappings
    /// and indices.
    pub fn new(
        annotations: GoAnnotations,
        ontology: &FullCsrOntology,
        pop_genes: Option<&HashSet<String>>,
    ) -> Self {
        // 1. Load and Filter Annotations (String-based)
        let raw_named_gene_to_term_sparse = get_annotation_map(&annotations);
        let named_gene_to_term_sparse =
            Self::filter_annotations(raw_named_gene_to_term_sparse, pop_genes);

        // 2. Build The Universe (Indices)
        let (gene_map, term_map) = Self::build_indices(&named_gene_to_term_sparse, ontology);

        // 3. Build all other Maps
        let mut named_term_to_gene_sparse: HashMap<TermId, HashSet<String>> = HashMap::new();
        let mut named_gene_to_term_dense: HashMap<String, HashSet<TermId>> = HashMap::new();
        let mut named_term_to_gene_dense: HashMap<TermId, HashSet<String>> = HashMap::new();

        for (gene_sym, terms) in &named_gene_to_term_sparse {
            // Ensure gene exists in dense map even if it has no terms
            named_gene_to_term_dense
                .entry(gene_sym.clone())
                .or_default(); // empty HashSet

            for term in terms {
                // A. Populate Sparse Reverse Map
                named_term_to_gene_sparse
                    .entry(term.clone())
                    .or_default()
                    .insert(gene_sym.clone());

                // B. Populate Dense Maps (True Path Rule)
                for ancestor in ontology.iter_term_and_ancestor_ids(term) {
                    named_term_to_gene_dense
                        .entry(ancestor.clone())
                        .or_default()
                        .insert(gene_sym.clone());

                    named_gene_to_term_dense
                        .entry(gene_sym.clone())
                        .or_default()
                        .insert(ancestor.clone());
                }
            }
        }

        // 5. Build Dense Matrices (Vec<Vec<usize>>) using the fixed indices
        // --- A. Term Matrices (Iterate term_map) ---
        let mut term_to_gene_dense = Vec::with_capacity(term_map.len());
        let mut term_to_gene_sparse = Vec::with_capacity(term_map.len());

        for term_id in &term_map {
            // Helper to convert Set<GeneSymbol> -> Sorted Vec<usize>
            let to_indices = |set: Option<&HashSet<String>>| -> IndexSet<usize> {
                match set {
                    Some(s) => {
                        let mut idxs: IndexSet<usize> = s
                            .iter()
                            .map(|g| gene_map.get_index_of(g).expect("Gene index missing"))
                            .collect();
                        idxs.sort_unstable();
                        idxs
                    }
                    None => IndexSet::new(), // Term exists in Dense but not Sparse -> Empty list
                }
            };

            term_to_gene_dense.push(to_indices(named_term_to_gene_dense.get(term_id)));
            term_to_gene_sparse.push(to_indices(named_term_to_gene_sparse.get(term_id)));
        }

        // --- B. Gene Matrices (Iterate gene_map) ---
        let mut gene_to_term_dense = Vec::with_capacity(gene_map.len());
        let mut gene_to_term_sparse = Vec::with_capacity(gene_map.len());

        for gene_sym in &gene_map {
            // Helper to convert Set<TermId> -> Sorted Vec<usize>
            let to_indices = |set: Option<&HashSet<TermId>>| -> IndexSet<usize> {
                match set {
                    Some(s) => {
                        let mut idxs: IndexSet<usize> = s
                            .iter()
                            .map(|t| term_map.get_index_of(t).expect("Term index missing"))
                            .collect();
                        idxs.sort_unstable();
                        idxs
                    }
                    None => IndexSet::new(),
                }
            };

            gene_to_term_dense.push(to_indices(named_gene_to_term_dense.get(gene_sym)));
            gene_to_term_sparse.push(to_indices(named_gene_to_term_sparse.get(gene_sym)));
        }

        Self {
            annotations,
            term_map,
            gene_map,
            gene_to_term_sparse,
            term_to_gene_sparse,
            gene_to_term_dense,
            term_to_gene_dense,
        }
    }

    /// Helper: Filter raw annotations based on population genes.
    fn filter_annotations(
        raw_map: HashMap<String, HashSet<TermId>>,
        pop_genes: Option<&HashSet<String>>,
    ) -> HashMap<String, HashSet<TermId>> {
        match pop_genes {
            Some(pop) => {
                let mut filtered_map = HashMap::new();
                for gene in pop {
                    // Only create an entry if the gene exists in the source annotations
                    if let Some(terms) = raw_map.get(gene) {
                        filtered_map.insert(gene.clone(), terms.clone());
                    }
                }
                filtered_map
            }
            None => raw_map,
        }
    }

    /// Helper: Collect all unique Genes and Terms (including inferred ancestors) to build IndexSets.
    fn build_indices(
        sparse_map: &HashMap<String, HashSet<TermId>>,
        ontology: &FullCsrOntology,
    ) -> (IndexSet<String>, IndexSet<TermId>) {
        let mut gene_map = IndexSet::new();
        let mut term_map = IndexSet::new();

        for (gene, terms) in sparse_map {
            gene_map.insert(gene.clone());
            for term in terms {
                term_map.insert(term.clone());
                // Pre-calculate ancestors so they have indices in the term_map
                for ancestor in ontology.iter_term_and_ancestor_ids(term) {
                    term_map.insert(ancestor.clone());
                }
            }
        }

        // Sorting ensures deterministic index assignment across runs
        gene_map.sort();
        term_map.sort();

        (gene_map, term_map)
    }

    // --- Accessors ---

    pub fn get_terms(&self) -> &IndexSet<TermId> {
        &self.term_map
    }

    pub fn get_term_by_index(&self, i: usize) -> &TermId {
        &self.term_map[i]
    }

    pub fn get_index_by_term(&self, term: &TermId) -> Option<usize> {
        self.term_map.get_index_of(term)
    }

    pub fn get_genes(&self) -> &IndexSet<String> {
        &self.gene_map
    }

    pub fn get_gene_by_index(&self, i: usize) -> &String {
        &self.gene_map[i]
    }

    pub fn get_index_by_gene(&self, gene: &String) -> Option<usize> {
        self.gene_map.get_index_of(gene)
    }

    pub fn get_gene_idxs_for_term_idx(&self, i: usize) -> &IndexSet<usize> {
        &self.term_to_gene_dense[i]
    }

    pub fn get_term_idxs_for_gene_idx(&self, i: usize) -> &IndexSet<usize> {
        &self.gene_to_term_dense[i]
    }

    pub fn get_terms_to_genes(&self, dense: bool) -> &Vec<IndexSet<usize>> {
        if dense {
            &self.term_to_gene_dense
        } else {
            &self.term_to_gene_sparse
        }
    }

    pub fn get_genes_to_terms(&self, dense: bool) -> &Vec<IndexSet<usize>> {
        if dense {
            &self.gene_to_term_dense
        } else {
            &self.gene_to_term_sparse
        }
    }
}
