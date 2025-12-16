use crate::bayesian::proposer::ToggleSwap;
use ToggleSwap::*;
use crate::core::AnnotationIndex;

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

pub(crate) struct MgsaState<'a> {
    // 1. Parameters / terms we want to infer
    terms: Vec<bool>, // maybe BitSet later
    n: usize,
    n_on: usize,
    n_off: usize,

    // 2. Latent space that connects parameters to observations
    latent : Vec<usize>,

    // 3. The Graph
    annotations: & 'a AnnotationIndex
}

impl<'a> MgsaState<'a> {
    pub(crate) fn new(terms: Vec<bool>, annotations: &'a AnnotationIndex) -> MgsaState {
        let n = terms.len();
        let n_on = terms.iter().filter(|&x| *x == true).count();
        let n_off = terms.iter().filter(|&x| *x == false).count();
        assert_eq!(n, annotations.terms().len(), "Initial terms vector size does not match AnnotationIndex");

        let mut latent = Self::construct_latent(&terms, annotations);

        MgsaState {
            terms,
            n,
            n_on,
            n_off,
            latent,
            annotations
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

    fn construct_latent(terms: &[bool], annotations: &AnnotationIndex) -> Vec<usize> {
        // 1. Initialize with correct size (using semicolon!)
        let mut latent = vec![0; annotations.genes().len()];

        // 2. Iterate and accumulate
        for (term_idx, &active) in terms.iter().enumerate() {
            if active {
                // Use the getter we fixed earlier
                let gene_idxs = annotations.get_gene_idxs_for_term_idx(term_idx);
                for &gene_i in gene_idxs {
                    latent[gene_i] += 1;
                }
            }
        }
        latent
    }

    fn update_latent(&mut self, term_idx: usize, enable: bool) {
        let gene_indices = self.annotations.get_gene_idxs_for_term_idx(term_idx);
        if enable{
            for &gene_i in gene_indices{
                self.latent[gene_i] += 1;
            }
        }
        else{
            for &gene_i in gene_indices{
                debug_assert!(self.latent[gene_i] > 0, "Latent count underflow!");
                self.latent[gene_i] -= 1;
            }
        }
    }
}

impl<'a> State for MgsaState<'a> {
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

impl<'a> CountableState for MgsaState<'a> {
    fn n_active(&self) -> usize {
        self.n_on
    }

    fn n_inactive(&self) -> usize {
        self.n_off
    }
}
