use oboannotation::go::{GoAnnotations, stats::get_annotation_map};
use ontolius::{
    TermId,
    ontology::{HierarchyWalks, csr::FullCsrOntology},
};

use indexmap::IndexSet;
use std::collections::{HashMap, HashSet};

// Contains all GO annotations and provides methods to build annotation maps specific to the study and population sets.
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
    gene_to_term_sparse: Vec<Vec<usize>>,
    term_to_gene_sparse: Vec<Vec<usize>>,
    // dense / including ancestors / explicit
    gene_to_term_dense: Vec<Vec<usize>>,
    term_to_gene_dense: Vec<Vec<usize>>,
}

impl AnnotationIndex {
    /// ```rust
    /// Creates a new instance of the data structure by processing gene-to-term annotations,
    /// applying optional population gene filtering, and constructing all relevant mappings
    /// and indices.
    ///
    /// # Parameters
    /// - `annotations`: An instance of `GoAnnotations` representing the initial annotations
    ///   that map genes to GO terms.
    /// - `ontology`: A reference to the `FullCsrOntology` object, which provides ontology information
    ///   such as term relationships and ancestor terms.
    /// - `pop_genes`: An optional reference to a `HashSet<String>` representing a population gene set.
    ///   If provided, only annotations for genes in this set will be included.
    ///
    /// # Returns
    /// An instance of the constructed data structure with precomputed mappings, indices, and matrices.
    ///
    /// # Steps
    /// 1. **Load Sparse Annotations**: Extracts a raw gene-to-term mapping from the provided annotations.
    /// 2. **Restrict to Population Gene Set**: If `pop_genes` is provided, filters the raw annotations
    ///    to include only genes present in the population gene set. If not provided, all genes are included.
    /// 3. **Build Sparse and Dense Maps**:
    ///    - `named_gene_to_term_sparse`: A sparse mapping from gene symbols to their directly annotated terms.
    ///    - `named_term_to_gene_sparse`: A reverse mapping (sparse) from terms to the genes annotated to them.
    ///    - `named_gene_to_term_dense`: A dense mapping that includes not only direct annotations but also
    ///      inferred annotations using the "True Path Rule" (including ancestors).
    ///    - `named_term_to_gene_dense`: A reverse dense mapping from terms (including inferred ancestors)
    ///      to genes.
    /// 4. **Build Index Maps**
    pub fn new(
        annotations: GoAnnotations,
        ontology: &FullCsrOntology,
        pop_genes: Option<&HashSet<String>>,
    ) -> Self {
        // 2. Load Sparse Annotations (implicit)
        let raw_gene_to_term_sparse = get_annotation_map(&annotations);

        // 3. Restrict to Population Gene Set
        let named_gene_to_term_sparse: HashMap<String, HashSet<TermId>> = match pop_genes {
            Some(pop) => {
                let mut filtered_map = HashMap::new();
                for gene in pop {
                    // Check if gene exists in annotations. If so, create an entry with
                    // key = GeneSymbol and value=Term list (cloned). If not, not create entry with
                    // an empty term list.
                    let terms = raw_gene_to_term_sparse
                        .get(gene)
                        .cloned()
                        .unwrap_or_default();
                    filtered_map.insert(gene.clone(), terms);
                }
                filtered_map
            }
            None => raw_gene_to_term_sparse,
        };

        // 4. Build all other Maps
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

        // 5. Build Index Maps. Sort for deterministic behavior across runs.
        let mut term_map: IndexSet<TermId> = named_term_to_gene_dense.keys().cloned().collect();
        term_map.sort(); // Deterministic ordering

        let mut gene_map: IndexSet<String> = named_gene_to_term_dense.keys().cloned().collect();
        gene_map.sort(); // Deterministic ordering

        // 6. Build Dense Matrices (Vec<Vec<usize>>) using the fixed indices
        // --- A. Term Matrices (Iterate term_map) ---
        let mut term_to_gene_dense = Vec::with_capacity(term_map.len());
        let mut term_to_gene_sparse = Vec::with_capacity(term_map.len());

        for term_id in &term_map {
            // Helper to convert Set<GeneSymbol> -> Sorted Vec<usize>
            let to_indices = |set: Option<&HashSet<String>>| -> Vec<usize> {
                match set {
                    Some(s) => {
                        let mut idxs: Vec<usize> = s
                            .iter()
                            .map(|g| gene_map.get_index_of(g).expect("Gene index missing"))
                            .collect();
                        idxs.sort_unstable();
                        idxs
                    }
                    None => Vec::new(), // Term exists in Dense but not Sparse -> Empty list
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
            let to_indices = |set: Option<&HashSet<TermId>>| -> Vec<usize> {
                match set {
                    Some(s) => {
                        let mut idxs: Vec<usize> = s
                            .iter()
                            .map(|t| term_map.get_index_of(t).expect("Term index missing"))
                            .collect();
                        idxs.sort_unstable();
                        idxs
                    }
                    None => Vec::new(),
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

    pub fn terms(&self) -> &IndexSet<TermId> {
        &self.term_map
    }

    /// Convert an index to its TermId
    pub fn get_index_term(&self, i: usize) -> &TermId {
        &self.term_map[i]
    }

    /// Convert a TermId to its index
    pub fn get_term_index(&self, term: &TermId) -> Option<usize> {
        self.term_map.get_index_of(term)
    }

    pub fn genes(&self) -> &IndexSet<String> {
        &self.gene_map
    }

    pub fn get_index_gene(&self, i: usize) -> &String {
        &self.gene_map[i]
    }

    pub fn get_gene_index(&self, gene: &String) -> Option<usize> {
        self.gene_map.get_index_of(gene)
    }

    pub fn get_gene_idxs_for_term_idx(&self, i: usize) -> &Vec<usize> {
        &self.term_to_gene_dense[i]
    }

    pub fn get_term_idxs_for_gene_idx(&self, i: usize) -> &Vec<usize> {
        &self.gene_to_term_dense[i]
    }

    pub fn get_terms_to_genes(&self, dense: bool) -> &Vec<Vec<usize>> {
        if dense {
            &self.term_to_gene_dense
        } else {
            &self.term_to_gene_sparse
        }
    }

    pub fn get_genes_to_terms(&self, dense: bool) -> &Vec<Vec<usize>> {
        if dense {
            &self.gene_to_term_dense
        } else {
            &self.gene_to_term_sparse
        }
    }
}
