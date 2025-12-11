use rand::Rng;
use crate::better_bayesian::state::{CountableState, State};
use ToggleSwap::{Swap, Toggle};

pub trait Proposer<S: State> {
    /// Generates a potential move
    fn propose<R: Rng>(&self, state: &S, rng: &mut R) -> S::Move;

    /// Calculates q(x'|x) / q(x|x')
    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> f64;
}


#[derive(Clone, Copy)]
pub enum ToggleSwap {
    Toggle(usize),
    Swap(usize, usize),
}

impl<S> Proposer<S> for ToggleSwap
where
    S : CountableState<Move = ToggleSwap, Value = bool>
{
    fn propose<R : Rng>(&self, state : &S, rng : &mut R) -> S::Move {
        // Every possible state transition is equally likely.
        let n = state.n_all();
        let na = state.n_active();
        let ni = state.n_inactive();

        let m = (n + na * ni);
        let x = rng.random_range(0..m);

        // In m cases we flip the state of a single term.
        if x < n {
            Toggle(x)
        }
        // in m_on * m_off cases we swap the states of two terms with different states.
        else {
            // map random number x to pairs of indices a, b
            let k = x - n;

            let (i, j) = find_swap_indices(k, state)
                .expect("Swap index out of range");

            Swap(i, j)
        }
    }

    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> f64 {
        match *m {
            Toggle(i) => {
                let n_current = (state.n_all() + state.n_active() * state.n_inactive()) as f64;

                // Calculate the change in possible moves (delta).

                // Active -> Inactive : new_n = current_n + n_on - n_off - 1
                // Inactive -> Active : new_n = current_n + n_off - n_on - 1
                let diff = state.n_active() as f64 - state.n_inactive() as f64;
                // If terms[i] is true (Active->Inactive), we use 'diff-1'.
                // If terms[i] is false (Inactive->Active), we use '-diff-1'.
                let delta = if state.get(i) { diff } else { -diff } - 1.0;

                let n_proposed = n_current + delta;
                (n_current / n_proposed).ln()
            }
            // Swapping preserves n_on and n_off, so the state space size N stays constant.
            // ln(N / N) = ln(1) = 0
            Swap(_, _) => 0.0,
        }
    }
}

/// Maps a random number `k` to indices (index_on, index_off) in the terms vector.
/// Assumption: k < m_on * m_off
fn find_swap_indices<S>(k : usize, state : &S) -> Option<(usize, usize)>
where
    S: State<Value = bool> + CountableState
{
    let n = state.n_all();
    let na = state.n_active();
    let ni = state.n_inactive();

    if k >= na * ni {
        return None;
    }

    let target_on_nth = k % na;
    let target_off_nth = k / na;

    let mut on_idx = None;
    let mut off_idx = None;
    let mut current_on_count = 0;
    let mut current_off_count = 0;

    // Scan the vector once to find the actual indices of "nth active" and "nth inactive" terms.
    for i in 0..n{
        if state.get(i) {
            if current_on_count == target_on_nth {
                on_idx = Some(i);
            }
            current_on_count += 1;
        } else {
            if current_off_count == target_off_nth {
                off_idx = Some(i);
            }
            current_off_count += 1;
        }

        if on_idx.is_some() && off_idx.is_some() {
            // FIX 2: Return immediately here for cleaner flow
            return Some((on_idx.unwrap(), off_idx.unwrap()));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::better_bayesian::state::Terms;
    #[test]
    fn test_swap_indices_mapping() {
        let terms = Terms::new(vec![true, false, false, true, false]);

        // Case k=0: target_on=0, target_off=0 -> Expect indices (0, 1)
        assert_eq!(find_swap_indices(0, &terms).unwrap(), (0, 1));

        // Case k=1: target_on=1, target_off=0 -> Expect indices (3, 1)
        assert_eq!(find_swap_indices(1, &terms).unwrap(), (3, 1));

        // Case k=2: target_on=0, target_off=1 -> Expect indices (0, 2)
        assert_eq!(find_swap_indices(2, &terms).unwrap(), (0, 2));

        // Case k=3: target_on=1, target_off=1 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(3,&terms).unwrap(), (3, 2));

        // Case k=4: target_on=0, target_off=2 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(4, &terms).unwrap(), (0, 4));

        // Case k=5: target_on=1, target_off=2 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(5, &terms).unwrap(), (3, 4));
    }

    #[test]
    #[should_panic]
    fn test_out_of_bounds() {
        // If we ask for a k that implies more active terms than exist
        let terms = Terms::new(vec![true, false]);
        let m_on = 1;
        let m_off = 1;

        find_swap_indices(5, &terms).unwrap();
    }
}