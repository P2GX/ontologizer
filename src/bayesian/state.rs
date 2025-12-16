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
pub trait CountableState: State<Value = bool> {
    fn n_active(&self) -> usize;
    fn n_inactive(&self) -> usize;
}

pub(crate) struct MgsaState {
    // 1. Parameters / terms we want to infer
    terms: Vec<bool>, // maybe BitSet later
    n: usize,
    n_on: usize,
    n_off: usize,

    // 2. Latent space between parameters to observations. Here it mimics observations.
    latent: Vec<usize>,

    // 3. The Network structure that connects terms to latent state, derived from Annotations
    terms_to_genes: Vec<Vec<usize>>,
}

impl MgsaState {
    pub(crate) fn new(
        terms: Vec<bool>,
        terms_to_genes: Vec<Vec<usize>>,
        n_genes: usize,
    ) -> MgsaState {
        let n = terms.len();
        let n_on = terms.iter().filter(|&x| *x == true).count();
        let n_off = terms.iter().filter(|&x| *x == false).count();
        assert_eq!(
            n,
            terms_to_genes.len(),
            "Initial terms vector size does not match network structure"
        );

        let mut latent = Self::construct_latent(&terms, &terms_to_genes, n_genes);

        MgsaState {
            terms,
            n,
            n_on,
            n_off,
            latent,
            terms_to_genes,
        }
    }

    /// Returns true if the cached counts match the actual vector data (debug only).
    fn check_consistency(&self) -> bool {
        let actual_on = self.terms.iter().filter(|&&t| t).count();
        let actual_off = self.terms.iter().filter(|&&t| !t).count();

        // Return result rather than panicking, so we can use it in assert!
        actual_on == self.n_on && actual_off == self.n_off
    }

    fn toggle_term(&mut self, term_idx: usize) {
        // Assess whether term_idx is activated or inactivated
        let is_becoming_active = !self.terms[term_idx];

        // Update terms
        self.terms[term_idx] = is_becoming_active;
        // Update n_on/n_off count
        if is_becoming_active {
            self.n_on += 1;
            self.n_off -= 1;
        } else {
            self.n_on -= 1;
            self.n_off += 1;
        }

        // Update Latent
        self.update_latent(term_idx, is_becoming_active);
    }

    fn construct_latent(
        terms: &[bool],
        terms_to_genes: &Vec<Vec<usize>>,
        n_genes: usize,
    ) -> Vec<usize> {
        // 1. Initialize with correct size
        let mut latent = vec![0; n_genes];

        // 2. Iterate and accumulate
        for (term_idx, &active) in terms.iter().enumerate() {
            if active {
                // Use the getter we fixed earlier
                let gene_indices = &terms_to_genes[term_idx];
                for &gene_i in gene_indices {
                    latent[gene_i] += 1;
                }
            }
        }
        latent
    }

    fn update_latent(&mut self, term_idx: usize, enable: bool) {
        let gene_indices = &self.terms_to_genes[term_idx];
        if enable {
            for &gene_i in gene_indices {
                self.latent[gene_i] += 1;
            }
        } else {
            for &gene_i in gene_indices {
                debug_assert!(self.latent[gene_i] > 0, "Latent count underflow!");
                self.latent[gene_i] -= 1;
            }
        }
    }
}

impl State for MgsaState {
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
            Toggle(i) => self.toggle_term(i),
            Swap(i, j) => {
                // A swap is just two toggles (On->Off, Off->On)
                self.toggle_term(i);
                self.toggle_term(j);
            }
        }
        debug_assert!(self.check_consistency());
    }

    /// Revert move and update n_on, n_off count
    fn revert(&mut self, m: &ToggleSwap) {
        self.apply(m)
    }
}

impl<'a> CountableState for MgsaState {
    fn n_active(&self) -> usize {
        self.n_on
    }

    fn n_inactive(&self) -> usize {
        self.n_off
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
    fn create_state(terms: Vec<bool>) -> MgsaState {
        // Hardcoded integer graph (no strings needed!)
        let terms_to_genes = vec![
            vec![0, 1], // T0 activates G0, G1
            vec![1, 2], // T1 activates G1, G2
        ];

        let n_genes = 3;

        // Verify the test writer didn't mess up the input vector size
        assert_eq!(
            terms.len(),
            terms_to_genes.len(),
            "Test config must have length 2"
        );

        MgsaState::new(terms, terms_to_genes, n_genes)
    }

    #[test]
    fn test_initialization_and_consistency() {
        // Initial state: T0 is ON, T1 is OFF
        let mut state = create_state(vec![true, false]);

        // 1. Check Counters
        assert_eq!(state.n_active(), 1);
        assert_eq!(state.n_inactive(), 1);
        assert!(state.check_consistency());

        // 2. Check Latent Space Calculation
        // T0 is ON -> Activates G0 and G1
        // T1 is OFF
        // Expected Latent: G0=1, G1=1, G2=0
        assert_eq!(state.latent, vec![1, 1, 0]);
    }

    #[test]
    fn test_toggle_updates() {
        let mut state = create_state(vec![true, false]);
        // --- Action: Toggle T1 ON ---
        // New State: T0=On, T1=On
        state.apply(&Toggle(1));

        assert_eq!(state.n_active(), 2);
        assert!(state.get(1)); // T1 should be true

        // Check Latent:
        // T0 adds: G0, G1
        // T1 adds: G1, G2
        // Result: G0=1, G1=2, G2=1
        assert_eq!(state.latent, vec![1, 2, 1]);
        assert!(state.check_consistency());

        // --- Action: Toggle T0 OFF ---
        // New State: T0=Off, T1=On
        state.apply(&Toggle(0));

        assert_eq!(state.n_active(), 1);
        assert!(!state.get(0)); // T0 should be false

        // Check Latent:
        // T0 removed: G0-1, G1-1
        // Result: G0=0, G1=1, G2=1
        assert_eq!(state.latent, vec![0, 1, 1]);
    }

    #[test]
    fn test_swap_move() {
        let mut state = create_state(vec![true, false]);

        // Swap(0, 1) should flip both
        // T0 -> Off, T1 -> On
        state.apply(&Swap(0, 1));

        assert!(!state.get(0));
        assert!(state.get(1));

        // Latent check:
        // Only T1 active (G1, G2)
        assert_eq!(state.latent, vec![0, 1, 1]);
    }

    #[test]
    fn test_revert() {
        let mut state = create_state(vec![true, false]);

        let m = Toggle(1);

        // Apply
        state.apply(&m);
        assert!(state.get(1));
        assert_eq!(state.latent, vec![1, 2, 1]);

        // Revert
        state.revert(&m);
        assert!(!state.get(1)); // Back to false
        assert_eq!(state.latent, vec![1, 1, 0]); // Back to original
        assert!(state.check_consistency());
    }
}
