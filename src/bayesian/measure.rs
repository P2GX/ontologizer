use crate::core::result::Measure;

#[derive(Debug)]
pub struct Probability {
    /// The probability score.
    probability: f64,
    /// Diagnostic metric (e.g., number of swaps or visits).
    swap: usize,
}

impl Probability {
    pub fn new(probability: f64, swap: usize) -> Self {
        Probability { probability, swap }
    }
}

impl Measure for Probability {
    fn score(&self) -> f64 {
        self.probability
    }

    fn diagnostics(&self) -> Option<String> {
        Some(self.swap.to_string())
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mean {
    /// The mean value.
    mean: f64,
    /// Diagnostic metric (e.g., number of swaps or visits).
    swap: usize,
}

impl Mean {
    pub fn new(mean: f64, swap: usize) -> Self {
        Mean { mean, swap }
    }
}

impl Measure for Mean {
    fn score(&self) -> f64 {
        self.mean
    }

    fn diagnostics(&self) -> Option<String> {
        Some(self.swap.to_string())
    }
}
