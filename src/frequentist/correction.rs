use serde::{Deserialize, Serialize};

pub trait Adjustment {
    fn adjust<T, F>(&self, items: &mut [T], extract: F)
    where
        F: Fn(&mut T) -> &mut f64;

    #[allow(dead_code)]
    fn name(&self) -> &'static str;
    // Used by the Ontologizer frontend
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(rename_all = "lowercase")]
pub enum Correction {
    Bonferroni,
    BonferroniHolm,
    BenjaminHochberg,
    None,
}

impl Correction {
    pub fn all() -> &'static [Correction] {
        &[
            Correction::Bonferroni,
            Correction::BonferroniHolm,
            Correction::BenjaminHochberg,
            Correction::None,
        ]
    }
}

impl Adjustment for Correction {
    fn adjust<T, F>(&self, items: &mut [T], extract: F)
    where
        F: Fn(&mut T) -> &mut f64,
    {
        match self {
            Correction::Bonferroni => {
                let n = items.len() as f64;
                for item in items.iter_mut() {
                    let p_ref = extract(item);
                    *p_ref = (*p_ref * n).min(1.0);
                }
            }
            Correction::BonferroniHolm => {
                let n = items.len();
                if n == 0 {
                    return;
                }

                let mut indices: Vec<usize> = (0..n).collect();
                let p_values_copy: Vec<f64> = items.iter_mut().map(|item| *extract(item)).collect();

                indices.sort_by(|&a, &b| {
                    p_values_copy[a]
                        .partial_cmp(&p_values_copy[b])
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                let m = n as f64;
                let mut max_prev = 0.0;

                for (rank_minus_1, &idx) in indices.iter().enumerate() {
                    let rank = rank_minus_1 as f64 + 1.0;
                    let correction_factor = m - rank + 1.0;
                    let p_ref = extract(&mut items[idx]);
                    let raw_p = *p_ref;
                    let mut adj_p = (raw_p * correction_factor).min(1.0);
                    adj_p = adj_p.max(max_prev);
                    max_prev = adj_p;
                    *p_ref = adj_p;
                }
            }
            Correction::BenjaminHochberg => {
                todo!("Benjamini-Hochberg correction not yet implemented")
            }
            Correction::None => {}
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Correction::Bonferroni => "Bonferroni",
            Correction::BonferroniHolm => "Bonferroni-Holm",
            Correction::BenjaminHochberg => "Benjamini-Hochberg",
            Correction::None => "No Correction",
        }
    }
}
