use crate::core::result::Measure;

pub struct PValue {
    pvalues: Vec<f64>,
    chose: Vec<(u32, u32, u32, u32)>,
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
        self.chose.iter().map(|&x| Some(format!("{:?}", x)))
    }

    fn get_score(&self, i: usize) -> f64 {
        self.pvalues[i]
    }

    fn get_diagnostics(&self, i: usize) -> Option<String> {
        Some(format!("{:?}", self.chose[i]))
    }
}
