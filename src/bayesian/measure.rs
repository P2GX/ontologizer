use crate::core::result::Measure;

#[derive(Debug)]
pub struct Probability {
    probabilities: Vec<f64>,
    swaps: Vec<usize>,
}

impl Probability {
    pub fn new(probabilities: Vec<f64>, swaps: Vec<usize>) -> Self {
        Probability {
            probabilities,
            swaps,
        }
    }
    pub fn iter(&self) -> std::slice::Iter<f64> {
        self.probabilities.iter()
    }
}

impl Measure for Probability {
    fn scores(&self) -> impl Iterator<Item = f64> {
        self.probabilities.iter().cloned()
    }

    fn diagnostics(&self) -> impl Iterator<Item = Option<String>> {
        self.swaps.iter().map(|&x| Some(x.to_string()))
    }

    fn get_score(&self, i: usize) -> f64 {
        self.probabilities[i] as f64
    }

    fn get_diagnostics(&self, i: usize) -> Option<f64> {
        Some(self.swaps[i] as f64)
    }
}
