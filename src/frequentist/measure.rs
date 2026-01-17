use crate::core::result::Measure;

pub struct PValue {
    pub pvalues : Vec<f64>,
}

impl PValue {
    pub fn iter(&self) -> std::slice::Iter<f64> {
        self.pvalues.iter()
    }
}

impl Measure for PValue {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.pvalues.iter().cloned()
    }

    fn diagnostics(&self) -> impl Iterator<Item = Option<String>> {
        std::iter::repeat(None)
    }

    fn get_score(&self, i: usize) -> f64 {
        self.pvalues[i] as f64
    }

    fn get_diagnostics(&self, i: usize) -> Option<f64> { None }
}