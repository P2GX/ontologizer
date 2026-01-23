///
pub trait Correction {
    fn adjust<T, F>(&self, items: &mut [T], extract: F)
    where
        F: Fn(&mut T) -> &mut f64;
    fn name(&self) -> &'static str;
}

pub struct Bonferroni;

impl Correction for Bonferroni {
    fn adjust<T, F>(&self, items: &mut [T], extract: F)
    where
        F: Fn(&mut T) -> &mut f64,
    {
        let n = items.len() as f64;

        // Iterate over items to access the p-values via the closure
        for item in items.iter_mut() {
            // Get mutable reference to the p-value
            let p_ref = extract(item);

            // Apply correction in-place
            *p_ref = (*p_ref * n).min(1.0);
        }
    }

    fn name(&self) -> &'static str {
        "Bonferroni"
    }
}

pub struct BonferroniHolm;

impl Correction for BonferroniHolm {
    fn adjust<T, F>(&self, items: &mut [T], extract: F)
    where
        F: Fn(&mut T) -> &mut f64,
    {
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

            // enforce monotony
            adj_p = adj_p.max(max_prev);
            max_prev = adj_p;

            // write back
            *p_ref = adj_p;
        }
    }

    fn name(&self) -> &'static str {
        "Bonferroni-Holm"
    }
}

pub struct None;

impl Correction for None {
    fn adjust<T, F>(&self, _pvalue: &mut [T], _extract: F)
    where
        F: Fn(&mut T) -> &mut f64,
    {
        return; // No adjustment is made
    }

    fn name(&self) -> &'static str {
        "No Correction"
    }
}
