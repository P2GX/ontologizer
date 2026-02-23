mod bayesian;
mod core;
mod frequentist;

pub use core::{AnnotationIndex, EnrichmentItem, EnrichmentResult, GeneSet, Ontology};

pub use bayesian::analysis::analysis as bayesian_analysis;

pub use frequentist::analysis::analysis as frequentist_analysis;
