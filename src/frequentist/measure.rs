use crate::core::result::Measure;

pub struct PValue {
    pub(crate) pvalue: f64,
    // (n_annotated_study (k), n_study (n), n_annotated_pop (K), n_pop (N))
    counts: (u32, u32, u32, u32),
}

impl PValue {
    pub fn new(pvalue: f64, counts: (u32, u32, u32, u32)) -> Self {
        Self { pvalue, counts }
    }
}

impl Measure for PValue {
    fn score(&self) -> f64 {
        self.pvalue
    }

    fn diagnostics(&self) -> Option<String> {
        Some(format!("{:?}", self.counts))
    }
}
