mod annotations;
mod geneset;
pub(crate) mod result;
mod problem;

pub use geneset::{GeneSet, load_gene_set, separate_gene_set};

pub use annotations::AnnotationIndex;
