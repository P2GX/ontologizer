use crate::better_bayesian::proposer::ToggleSwap;
use ToggleSwap::*;

// A trait that guarantees that *STATE* knows how to sample itself by drawing and applying *MOVE*.
pub trait State {
    type Move;
    type Value: Copy;
    fn get(&self, i: usize) -> Self::Value;
    fn n_all(&self) -> usize;
    fn apply(&mut self, m: &Self::Move);
    fn revert(&mut self, m: &Self::Move);
}

// A trait that guarantees that *STATE* takes boolean values (activer/inactive).
pub trait CountableState: State<Value = bool> {
    fn n_active(&self) -> usize;
    fn n_inactive(&self) -> usize;
}

pub(crate) struct Terms {
    // Computational vector used by *ALGORITHM*. Mapping term indices to names by *ANNOTATIONS**
    terms: Vec<bool>, // maybe BitSet later
    n: usize,
    n_on: usize,
    n_off: usize,
}

impl Terms {
    pub(crate) fn new(terms: Vec<bool>) -> Terms {
        let n = terms.len();
        let n_on = terms.iter().filter(|&x| *x == true).count();
        let n_off = terms.iter().filter(|&x| *x == false).count();
        Terms {
            terms,
            n,
            n_on,
            n_off,
        }
    }

    /// Returns true if the cached counts match the actual vector data (debug only).
    fn check_consistency(&self) -> bool {
        let actual_on = self.terms.iter().filter(|&&t| t).count();
        let actual_off = self.terms.iter().filter(|&&t| !t).count();

        // Return result rather than panicking, so we can use it in assert!
        actual_on == self.n_on && actual_off == self.n_off
    }
}

impl State for Terms {
    type Move = ToggleSwap;

    type Value = bool;

    fn get(&self, i: usize) -> bool {
        self.terms[i]
    }

    fn n_all(&self) -> usize {
        self.n
    }
    /// Revert move and update n_on, n_off count
    fn apply(&mut self, m: &ToggleSwap) {
        match *m {
            Toggle(i) => {
                let is_now_on = !self.terms[i];
                if is_now_on {
                    self.n_on += 1;
                    self.n_off -= 1;
                } else {
                    self.n_on -= 1;
                    self.n_off += 1;
                }

                self.terms[i] = is_now_on
            }
            Swap(i_on, i_off) => {
                self.terms[i_on] = false;
                self.terms[i_off] = true;
            }
        }
        debug_assert!(self.check_consistency());
    }

    /// Revert move and update n_on, n_off count
    fn revert(&mut self, m: &ToggleSwap) {
        match *m {
            Toggle(i) => {
                let is_now_on = !self.terms[i];
                if is_now_on {
                    self.n_on += 1;
                    self.n_off -= 1;
                } else {
                    self.n_on -= 1;
                    self.n_off += 1;
                }

                self.terms[i] = is_now_on
            }
            Swap(i_on, i_off) => {
                self.terms[i_on] = true;
                self.terms[i_off] = false;
            }
        }
        debug_assert!(self.check_consistency());
    }
}

impl CountableState for Terms {
    fn n_active(&self) -> usize {
        self.n_on
    }

    fn n_inactive(&self) -> usize {
        self.n_off
    }
}
