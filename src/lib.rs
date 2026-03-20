mod bayesian;
mod core;
mod frequentist;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "method")]
pub enum Method {
    Frequentist {
        background: Background,
        correction: Correction,
    },
    Bayesian,
}

pub use frequentist::{Background, Correction};

pub use core::{AnalysisResult, AnnotationIndex, EnrichmentItem, GeneSet, Ontology};

pub use bayesian::analysis::analysis as bayesian_analysis;

pub use frequentist::analysis::analysis as frequentist_analysis;
