/// Represents a discrete change in the Term activation state.
#[derive(Clone, Copy, Debug)]
pub enum ToggleSwap {
    Toggle(usize),
    Swap(usize, usize),
}

/// Represents a continuous change in a parameter value.
#[derive(Clone, Copy, Debug)]
pub struct Increment {
    pub index: usize,
    pub delta: f64,
}

/// A union of possible moves in the MGSA state space.
#[derive(Clone, Debug)]
pub enum MgsaMove {
    Term(ToggleSwap),
    Parameter(Increment),
}

/// A trait defining the minimal requirements for a state in the Metropolis-Hastings algorithm.
///
/// A state must be able to apply a move to transition to a new state,
/// and revert that move to return to the previous state.
pub trait State {
    type Move;
    fn apply(&mut self, m: &Self::Move);
    fn revert(&mut self, m: &Self::Move);
}

// ==========================================
// TERM STATE
// ==========================================

/// Manages the configuration of active terms.
///
/// It maintains efficient indices to allow O(1) sampling of active/inactive terms.
pub(crate) struct TermState {
    terms: Vec<bool>,             // Boolean vector indicating if term `i` is active.
    active_indices: Vec<usize>,   // List of indices where `terms[i]` is true.
    inactive_indices: Vec<usize>, // List of indices where `terms[i]` is false.
    term_map: Vec<usize>,
    // A map pointing to the position of each term index in either `active_indices` or `inactive_indices`.
    // `active_indices[term_map[i]] == i` if `terms[i]` is true.
}

impl TermState {
    pub(crate) fn new(terms: Vec<bool>) -> TermState {
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

        TermState {
            terms,
            active_indices,
            inactive_indices,
            term_map,
        }
    }

    pub(crate) fn n_terms(&self) -> usize {
        self.terms.len()
    }

    pub fn get(&self, i: usize) -> bool {
        self.terms[i]
    }

    pub(crate) fn n_active(&self) -> usize {
        self.active_indices.len()
    }
    pub(crate) fn n_inactive(&self) -> usize {
        self.inactive_indices.len()
    }
    pub(crate) fn get_active(&self, k: usize) -> usize {
        self.active_indices[k]
    }
    pub(crate) fn get_inactive(&self, k: usize) -> usize {
        self.inactive_indices[k]
    }

    /// Toggles the state of a term (On <-> Off) and updates internal indices
    fn toggle(&mut self, term_idx: usize) {
        let is_becoming_active = !self.terms[term_idx];
        self.terms[term_idx] = is_becoming_active;

        if is_becoming_active {
            // Move from Inactive list to Active list
            self.move_index(term_idx, false);
        } else {
            // Move from Active list to Inactive list
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

        let source_pos = self.term_map[term_idx]; // index in source that has to be removed

        // `remove` removes one element and shifts over the remaining elements, O(n)
        // `swap_remove` moves the last element to the position of the removed element, O(1)
        // The drawback: we must update the map for that moved element.
        let last_element = source[source.len() - 1];
        source.swap_remove(source_pos);

        // Update map for the element that was swapped into the gap (unless we removed the last one)
        if source_pos < source.len() {
            self.term_map[last_element] = source_pos;
        }

        // Add to Destination
        self.term_map[term_idx] = dest.len();
        dest.push(term_idx);
    }

    #[cfg(test)]
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
}

impl State for TermState {
    type Move = ToggleSwap;

    fn apply(&mut self, m: &ToggleSwap) {
        match *m {
            ToggleSwap::Toggle(i) => self.toggle(i),
            ToggleSwap::Swap(i, j) => {
                // A swap is just two toggles (On->Off, Off->On)
                self.toggle(i);
                self.toggle(j);
            }
        }
    }

    fn revert(&mut self, m: &ToggleSwap) {
        // Reverting a toggle/swap is the same as applying it again
        self.apply(m)
    }
}

// ==========================================
// PARAMETER STATE
// ==========================================
#[derive(Clone, Debug)]
pub struct ParameterState {
    p: f64,
    alpha: f64,
    beta: f64,
}

impl ParameterState {
    pub fn new(p: f64, alpha: f64, beta: f64) -> Self {
        Self { p, alpha, beta }
    }

    pub fn n_params(&self) -> usize {
        3
    }

    pub fn p(&self) -> f64 {
        self.p
    }
    pub fn alpha(&self) -> f64 {
        self.alpha
    }
    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn get(&self, index: usize) -> f64 {
        match index {
            0 => self.p,
            1 => self.alpha,
            2 => self.beta,
            _ => panic!("Invalid parameter index: {}", index),
        }
    }

    pub fn update(&mut self, index: usize, delta: f64) {
        match index {
            0 => self.p += delta,
            1 => self.alpha += delta,
            2 => self.beta += delta,
            _ => panic!("Invalid parameter index: {}", index),
        }
    }
}

impl State for ParameterState {
    type Move = Increment;

    fn apply(&mut self, m: &Self::Move) {
        self.update(m.index, m.delta);
    }

    fn revert(&mut self, m: &Self::Move) {
        self.update(m.index, -m.delta);
    }
}

// ==========================================
// MGSA COMPOSITE STATE
// ==========================================

/// The full state of the MGSA algorithm, combining discrete Terms and continuous Parameters.
pub struct MgsaState {
    pub terms: TermState,
    pub params: ParameterState,
}

impl MgsaState {
    pub fn new(terms: Vec<bool>, p: f64, alpha: f64, beta: f64) -> Self {
        Self {
            terms: TermState::new(terms),
            params: ParameterState::new(p, alpha, beta),
        }
    }

    pub(crate) fn n_terms(&self) -> usize {
        self.terms.n_terms()
    }

    pub(crate) fn n_terms_active(&self) -> usize {
        self.terms.active_indices.len()
    }
    pub(crate) fn n_terms_inactive(&self) -> usize {
        self.terms.inactive_indices.len()
    }
}

impl State for MgsaState {
    type Move = MgsaMove;

    fn apply(&mut self, m: &Self::Move) {
        match m {
            MgsaMove::Term(ts) => self.terms.apply(ts),
            MgsaMove::Parameter(inc) => self.params.apply(inc),
        }
    }

    fn revert(&mut self, m: &Self::Move) {
        match m {
            MgsaMove::Term(ts) => self.terms.revert(ts),
            MgsaMove::Parameter(inc) => self.params.revert(inc),
        }
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
            let state = TermState::new(vec![true, false]);

            assert_eq!(state.n_active(), 1);
            assert_eq!(state.n_inactive(), 1);

            // Verify indices
            assert_eq!(state.get_active(0), 0); // Index 0 is active
            assert_eq!(state.get_inactive(0), 1); // Index 1 is inactive

            assert!(state.check_consistency());
        }

        #[test]
        fn test_toggle_updates() {
            let mut state = TermState::new(vec![true, false]);

            // --- Action: Toggle T1 ON ---
            state.apply(&ToggleSwap::Toggle(1));

            assert_eq!(state.n_active(), 2);
            assert!(state.get(1)); // T1 should be true
            assert!(state.check_consistency());

            // --- Action: Toggle T0 OFF ---
            state.apply(&ToggleSwap::Toggle(0));

            assert_eq!(state.n_active(), 1);
            assert!(!state.get(0)); // T0 should be false
            assert!(state.check_consistency());
        }

        #[test]
        fn test_swap_move() {
            let mut state = TermState::new(vec![true, false]);

            // Swap(0, 1) should flip both: T0->Off, T1->On
            state.apply(&ToggleSwap::Swap(0, 1));

            assert!(!state.get(0));
            assert!(state.get(1));
            assert!(state.check_consistency());
        }
    }
}
