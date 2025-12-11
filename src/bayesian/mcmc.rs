use ontolius::TermId;
use rand::{Rng, rngs::ThreadRng};
use std::collections::HashMap;

pub enum TermMove {
    Flip(TermId),
    Swap(TermId, TermId),
}

/// A sampler for the Bayesian network. Proposes a new state based on the current state.
pub(crate) struct Sampler {
    rng: ThreadRng,
}

impl Sampler {
    pub(crate) fn new() -> Self {
        // todo!(should allow for a seed)
        Self { rng: rand::rng() }
    }

    /// Takes a Term configuartion T = {T1 = 0, T2 = 1, T3 = 1, ... Tm = 1} draws a move.
    /// A move is either a flip of a single term Ti or a swap of two terms Ti, Tj with Ti != Tj.
    pub fn draw_move(&mut self, terms: &HashMap<TermId, bool>) -> TermMove {
        let m = terms.len();
        let m0 = terms.iter().filter(|&(k, &v)| v == false).count();
        let m1 = terms.iter().filter(|&(k, &v)| v == true).count();

        // Every possible state transition is equally likely.
        let n = m + m0 * m1;
        let x = self.rng.random_range(0..n);

        // In m cases we flip the state of a single term.
        if x < m {
            let tid = terms.keys().nth(x).unwrap();
            TermMove::Flip(tid.clone())
        }
        // in m1*m0 cases we exchange the states of two terms with different states.
        else {
            // map random number x to pairs of indices a, b.
            let k = x - m;
            let a = k / m0; // 0..m1-1
            let b = k % m0; // 0..m0-1

            let mut on_left = a;
            let mut off_left = b;
            let mut on_pick = None;
            let mut off_pick = None;

            for (tid, &on) in terms.iter() {
                if on {
                    if on_left == 0 {
                        on_pick = Some(tid.clone());
                    } else {
                        on_left -= 1;
                    }
                } else {
                    if off_left == 0 {
                        off_pick = Some(tid.clone());
                    } else {
                        off_left -= 1;
                    }
                }
                if on_pick.is_some() && off_pick.is_some() {
                    break;
                }
            }
            TermMove::Swap(on_pick.unwrap(), off_pick.unwrap())
        }
    }

    /// Apply a move in place.
    pub fn apply_move(&self, terms: &mut HashMap<TermId, bool>, mv: &TermMove) {
        match mv {
            TermMove::Flip(t) => {
                if let Some(v) = terms.get_mut(&t) {
                    *v = !*v;
                }
            }
            TermMove::Swap(t_on, t_off) => {
                if let Some(v_on) = terms.get_mut(&t_on) {
                    *v_on = false;
                }
                if let Some(v_off) = terms.get_mut(&t_off) {
                    *v_off = true;
                }
            }
        }
    }

    /// Revert a previously applied move.
    pub fn revert_move(&self, terms: &mut HashMap<TermId, bool>, mv: &TermMove) {
        match mv {
            TermMove::Flip(t) => {
                if let Some(v) = terms.get_mut(&t) {
                    *v = !*v;
                }
            }
            TermMove::Swap(t_on, t_off) => {
                if let Some(v_on) = terms.get_mut(&t_on) {
                    *v_on = true;
                }
                if let Some(v_off) = terms.get_mut(&t_off) {
                    *v_off = false;
                }
            }
        }
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
