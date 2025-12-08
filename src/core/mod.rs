mod annotations;
mod geneset;
mod ontology;
mod problem;
mod result;

pub use geneset::{GeneSet, GeneSymbol, load_gene_set, separate_gene_set};

pub use annotations::AnnotationIndex;

pub use ontology::Ontologizer;
