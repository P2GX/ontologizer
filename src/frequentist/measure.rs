use crate::core::result::Measure;

pub struct PValue {
    pvalue: f64,
    chose: (u32, u32, u32, u32),
}

impl Measure for PValue {
    fn score(&self) -> f64 {
        self.pvalue
    }

    fn diagnostics(&self) -> Option<String> {
        Some(format!("{:?}", self.chose))
    }
}
