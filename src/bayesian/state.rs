use crate::bayesian::proposer::ToggleSwap;
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

// A trait that guarantees that *STATE* takes boolean values (active/inactive).
pub trait BinaryParameterState: State<Value = bool> {
    fn get_active(&self, k: usize) -> usize;
    fn get_inactive(&self, k: usize) -> usize;
    fn n_active(&self) -> usize;
    fn n_inactive(&self) -> usize;
}

pub(crate) struct MgsaState {
    // Parameters
    terms: Vec<bool>, // maybe BitSet later

    active_indices: Vec<usize>,
    inactive_indices: Vec<usize>,

    term_map: Vec<usize>,
}

impl MgsaState {
    pub(crate) fn new(terms: Vec<bool>) -> MgsaState {
        let n = terms.len();

        let mut active_indices = Vec::with_capacity(n);
        let mut inactive_indices = Vec::with_capacity(n);
        let mut term_map = vec![0; n];

        for (i, &is_active) in terms.iter().enumerate() {
            if is_active {
                term_map[i] = active_indices.len();
                active_indices.push(i);
            } else {
                term_map[i] = inactive_indices.len();
                inactive_indices.push(i);
            }
        }

        MgsaState {
            terms,
            active_indices,
            inactive_indices,
            term_map,
        }
    }

    /// Returns true if the cached counts and lists match the actual vector data (debug only).
    fn check_consistency(&self) -> bool {
        let n_on_vec = self.terms.iter().filter(|&&t| t).count();
        let n_off_vec = self.terms.iter().filter(|&&t| !t).count();

        // Check counts
        if n_on_vec != self.active_indices.len() || n_off_vec != self.inactive_indices.len() {
            return false;
        }

        // Check Map Integrity
        for (idx, &term_idx) in self.active_indices.iter().enumerate() {
            if !self.terms[term_idx] || self.term_map[term_idx] != idx {
                return false;
            }
        }
        for (idx, &term_idx) in self.inactive_indices.iter().enumerate() {
            if self.terms[term_idx] || self.term_map[term_idx] != idx {
                return false;
            }
        }

        true
    }

    // Swap is implemented as two toggles.
    fn toggle_term(&mut self, term_idx: usize) {
        let is_becoming_active = !self.terms[term_idx];
        self.terms[term_idx] = is_becoming_active;

        // Update of Active/Inactive lists
        if is_becoming_active {
            // Move from Inactive -> Active
            self.move_index(term_idx, false);
        } else {
            // Move from Active -> Inactive
            self.move_index(term_idx, true);
        }
    }

    /// A loco-ly optimized function to move index from one vec to another. Credits to Chatty.
    fn move_index(&mut self, term_idx: usize, from_active: bool) {
        let (source, dest) = if from_active {
            (&mut self.active_indices, &mut self.inactive_indices)
        } else {
            (&mut self.inactive_indices, &mut self.active_indices)
        };
        // We are using four vectors here:
        // 1. `terms`: the value of each term, [true, true, true, false, false]
        // 2. `active_indices`: the indices of *active* terms in `terms`, [0, 1, 2]
        // 3. `inactive_indices`: the indices of *inactive* terms in `terms`, [3, 4]
        // 4. `term_map`: the position of each term index in active / inactive, [0, 1, 2, 0, 1]
        // Let last map is crucial because it omits the need to scan active_indices/inactive_indices
        // to find the position of a certain term_idx i and allows access directly by i=term_map[k].

        let source_idx = self.term_map[term_idx]; // index in source that has to be removed

        // `remove` removes one element and shifts over the remaining elements, O(n)
        // `swap_remove` moves the last element to the position of the removed element, O(1)
        // The drawback: we must update the map for that moved element.
        let last_element = source[source.len() - 1];
        source.swap_remove(source_idx);

        // Update map for the element that was swapped into the gap (unless we removed the last one)
        if source_idx < source.len() {
            self.term_map[last_element] = source_idx;
        }

        // Add to Destination
        self.term_map[term_idx] = dest.len();
        dest.push(term_idx);
    }
}

impl State for MgsaState {
    type Move = ToggleSwap;
    type Value = bool;

    fn get(&self, i: usize) -> bool {
        self.terms[i]
    }

    fn n_all(&self) -> usize {
        self.terms.len()
    }
    /// Revert move and update n_on, n_off count
    fn apply(&mut self, m: &ToggleSwap) {
        match *m {
            Toggle(i) => self.toggle_term(i),
            Swap(i, j) => {
                // A swap is just two toggles (On->Off, Off->On)
                self.toggle_term(i);
                self.toggle_term(j);
            }
        }
        // debug_assert!(self.check_consistency());
    }

    /// Revert move and update n_on, n_off count
    fn revert(&mut self, m: &ToggleSwap) {
        self.apply(m)
    }
}

impl BinaryParameterState for MgsaState {
    fn get_active(&self, k: usize) -> usize {
        self.active_indices[k]
    }

    fn get_inactive(&self, k: usize) -> usize {
        self.inactive_indices[k]
    }
    fn n_active(&self) -> usize {
        self.active_indices.len()
    }

    fn n_inactive(&self) -> usize {
        self.inactive_indices.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to construct a dummy AnnotationIndex in memory
    /// We simulate:
    /// 3 Genes: G0, G1, G2
    /// 2 Terms: T0, T1
    /// Relations:
    ///   T0 -> [G0, G1]
    ///   T1 -> [G1, G2]
    ///
    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_initialization_and_consistency() {
            // T0=On, T1=Off
            let state = MgsaState::new(vec![true, false]);

            assert_eq!(state.n_active(), 1);
            assert_eq!(state.n_inactive(), 1);

            // Verify indices
            assert_eq!(state.get_active(0), 0); // Index 0 is active
            assert_eq!(state.get_inactive(0), 1); // Index 1 is inactive

            assert!(state.check_consistency());
        }

        #[test]
        fn test_toggle_updates() {
            let mut state = MgsaState::new(vec![true, false]);

            // --- Action: Toggle T1 ON ---
            state.apply(&Toggle(1));

            assert_eq!(state.n_active(), 2);
            assert!(state.get(1)); // T1 should be true
            assert!(state.check_consistency());

            // --- Action: Toggle T0 OFF ---
            state.apply(&Toggle(0));

            assert_eq!(state.n_active(), 1);
            assert!(!state.get(0)); // T0 should be false
            assert!(state.check_consistency());
        }

        #[test]
        fn test_swap_move() {
            let mut state = MgsaState::new(vec![true, false]);

            // Swap(0, 1) should flip both: T0->Off, T1->On
            state.apply(&Swap(0, 1));

            assert!(!state.get(0));
            assert!(state.get(1));
            assert!(state.check_consistency());
        }
    }
}
