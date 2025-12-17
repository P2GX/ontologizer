use crate::bayesian::state::{BinaryParameterState, State};
use ToggleSwap::{Swap, Toggle};
use rand::Rng;

/// Allows to sample from the sample space and to find the corresponding transition probability.
pub trait Proposer<S: State> {
    /// Generates move x -> x'
    fn propose<R: Rng>(&self, state: &S, rng: &mut R) -> S::Move;

    /// Calculates q(x'|x) / q(x|x')
    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> f64;
}

/// Explore the sample space by flipping the value of one element,
/// or by exchanging the (different) values of two elements.
#[derive(Clone, Copy)]
pub enum ToggleSwap {
    Toggle(usize),
    Swap(usize, usize),
}

pub struct UniformProposer;
impl UniformProposer {
    pub fn new() -> Self {
        Self
    }
}

impl<S> Proposer<S> for UniformProposer
where
    S: BinaryParameterState<Move = ToggleSwap, Value = bool>,
{
    /// Draw a random number representing the index to all possible toggles and swaps and propose
    /// the corresponding toggle or swap.
    fn propose<R: Rng>(&self, state: &S, rng: &mut R) -> S::Move {
        // Every possible state transition is equally likely.
        let n = state.n_all();
        let na = state.n_active();
        let ni = state.n_inactive();

        let m = n + na * ni;
        let x = rng.random_range(0..m);

        // In m cases we flip the state of a single term.
        if x < n {
            Toggle(x)
        }
        // in m_on * m_off cases we swap the states of two terms with different states.
        else {
            // map random number x to pairs of indices a, b
            let k = x - n;
            let i = state.get_kth_active(k % na);
            let j = state.get_kth_inactive(k / na);
            Swap(i, j)
        }
    }

    ///
    fn log_proposal_ratio(&self, state: &S, m: &S::Move) -> f64 {
        match *m {
            Toggle(i) => {
                let n_current = (state.n_all() + state.n_active() * state.n_inactive()) as f64;

                // Calculate the change in possible moves (delta).

                // Active -> Inactive: new_n = current_n + n_on - n_off - 1
                // Inactive -> Active: new_n = current_n + n_off - n_on - 1
                let diff = state.n_active() as f64 - state.n_inactive() as f64;
                // If terms[i] is true (Active->Inactive), we use 'diff-1'.
                // If terms[i] is false (Inactive->Active), we use '-diff-1'.
                let delta = if state.get(i) { diff } else { -diff } - 1.0;

                let n_proposed = n_current + delta;
                (n_proposed / n_current).ln()
            }
            // Swapping preserves n_on and n_off, so the state space size N stays constant.
            // ln(N / N) = ln(1) = 0
            Swap(_, _) => 0.0,
        }
    }
}
