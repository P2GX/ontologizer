mod bayesian;
mod core;
mod frequentist;

pub use core::{AnnotationIndex, GeneSet};

pub use bayesian::run::analysis as bayesian_analysis;

pub use frequentist::run::analysis as frequentist_analysis;
