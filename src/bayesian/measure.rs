use crate::core::result::Measure;

#[derive(Debug)]
pub struct Probability {
    probability: f64,
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
