mod annotations;
mod genes;
mod ontology;
pub(crate) mod result;

pub use ontology::Ontology;

pub use annotations::AnnotationIndex;

pub use genes::GeneSet;

pub use result::{EnrichmentItem, EnrichmentResult};
