use crate::AnnotationIndex;
use indexmap::IndexSet;
use ontolius::ontology::HierarchyWalks;
use ontolius::ontology::csr::FullCsrOntology;
use serde::{Deserialize, Serialize};

pub trait Restrict {
    fn restrict(
        &self,
        term_idx: usize,
        annotations: &AnnotationIndex,
        ontology: &FullCsrOntology,
    ) -> IndexSet<usize>;

    #[allow(dead_code)]
    fn name(&self) -> &'static str;
    // Used by the Ontologizer frontend
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "lowercase")]
pub enum Topology {
    Standard,
    ParentUnion,
    ParentIntersection,
}

impl Topology {
    pub fn all() -> &'static [Topology] {
        &[
            Topology::Standard,
            Topology::ParentUnion,
            Topology::ParentIntersection,
        ]
    }
}

impl Restrict for Topology {
    fn restrict(
        &self,
        term_idx: usize,
        annotations: &AnnotationIndex,
        ontology: &FullCsrOntology,
    ) -> IndexSet<usize> {
        match self {
            Topology::Standard => (0..annotations.get_genes().len()).collect(),
            Topology::ParentUnion => {
                let term_id = annotations.get_term_by_idx(term_idx);

                let parent_gene_sets: Vec<&IndexSet<usize>> = ontology
                    .iter_parent_ids(term_id)
                    .filter_map(|parent_id| annotations.get_idx_by_term(&parent_id))
                    .map(|parent_idx| annotations.get_gene_idxs_for_term_idx(parent_idx))
                    .collect();

                if parent_gene_sets.is_empty() {
                    return (0..annotations.get_genes().len()).collect();
                }

                let mut union = parent_gene_sets[0].clone();
                for set in &parent_gene_sets[1..] {
                    union.extend(set.iter().copied());
                }
                union
            }
            Topology::ParentIntersection => {
                let term_id = annotations.get_term_by_idx(term_idx);
                let parent_gene_sets: Vec<&IndexSet<usize>> = ontology
                    .iter_parent_ids(term_id)
                    .filter_map(|parent_id| annotations.get_idx_by_term(&parent_id))
                    .map(|parent_idx| annotations.get_gene_idxs_for_term_idx(parent_idx))
                    .collect();

                if parent_gene_sets.is_empty() {
                    return IndexSet::new();
                }

                let mut intersection = parent_gene_sets[0].clone();
                for set in &parent_gene_sets[1..] {
                    intersection.retain(|gene| set.contains(gene));
                }
                intersection
            }
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Topology::Standard => "Standard",
            Topology::ParentUnion => "Parent Union",
            Topology::ParentIntersection => "Parent Intersection",
        }
    }
}
