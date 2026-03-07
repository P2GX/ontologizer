mod bayesian;
mod core;
mod frequentist;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "method", rename_all = "lowercase")]
pub enum Method {
    Frequentist {
        topology: Topology,
        correction: Correction,
    },
    Bayesian,
}

pub use frequentist::{Correction, Topology};

pub use core::{AnalysisResult, AnnotationIndex, EnrichmentItem, GeneSet, Ontology};

pub use bayesian::analysis::analysis as bayesian_analysis;

pub use frequentist::analysis::analysis as frequentist_analysis;
