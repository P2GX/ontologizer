use indexmap::IndexMap;
use rand::{rng, Rng, rngs::ThreadRng};
use ontolius::TermId;


/// A sampler for the Bayesian network. Proposes a new state based on the current state.
struct Sampler{
    rng : ThreadRng
}

impl Sampler{
    fn new() -> Self{
        Self { rng : rng() }
    }

    /// Takes a Term configuartion T = {T1 = 0, T2 = 1, T3 = 1, ... Tm = 1} and
    /// samples a new configuration by flipping the state of one term Ti or
    /// exchanging the states of two terms Ti, Tj with different states Ti != Tj.
    pub fn sample(&mut self, current: &IndexMap<TermId, bool>) -> IndexMap<TermId, bool> {
        let mut state = current.clone();

        let m = state.len();
        let m0 = state.iter().filter(|&(k, &v)| v == false).count();
        let m1 = state.iter().filter(|&(k, &v)| v == true).count();

        // Every possible state transition is equally likely.
        let xi = self.rng.random_range(0..m + m0*m1);

        // In m cases we flip the state of a single term.
        if xi < m {
            let i = xi;
            // Flip the state of term i.
            if let Some((_, v)) = state.get_index_mut(i) {
                *v = !*v;
            }
        }
        // in m1*m0 cases we exchange the states of two terms with different states.
        else {
            // map random number xi to pairs of indices a, b.
            let i = (xi - m) / m0;
            let j = (xi - m) % m0;
            let mut on_idx = Vec::new();
            let mut off_idx = Vec::new();
            for (i, (_, &v)) in state.iter().enumerate() {
                if v { on_idx.push(i) } else { off_idx.push(i) }
            }
            let a = on_idx[i];
            let b = off_idx[j];

            // swap the states of a and b.
            if let Some((_, v)) = state.get_index_mut(a) {
                *v = !*v;
            }
            if let Some((_, v)) = state.get_index_mut(b) {
                *v = !*v;
            }
        }

        state
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sampler() {
        unimplemented!("The process, you must trust.")
    }

}