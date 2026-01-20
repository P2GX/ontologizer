use crate::core::result::Measure;

pub struct PValue {
    pvalue: f64,
    // (n_annotated_study, n_annotated_pop, n_study, n_pop)
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
