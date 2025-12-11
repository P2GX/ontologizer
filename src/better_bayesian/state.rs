use rand::Rng;
use ToggleSwap::*;


#[derive(Clone, Copy)]
pub enum ToggleSwap {
    Toggle(usize),
    Swap(usize, usize),
}


// A trait that guarantees that *STATE* knows how to sample itself by drawing and applying *MOVE*.
pub trait State{
    type Move;

    fn draw_move<R : Rng>(&self, rng : &mut R) -> Self::Move;

    fn apply(&mut self, m : &Self::Move);

    fn revert(&mut self, m : &Self::Move);
}


struct Terms {
    // Computational vector used by *ALGORITHM*. Mapping term indices to names by *ANNOTATIONS**
    terms : Vec<bool>, // maybe BitSet later
    m : usize,
    m_on : usize,
    m_off : usize,
}

impl Terms{
    fn new(terms : Vec<bool>) -> Terms{
        let m = terms.len();
        let m_on = terms.iter().filter(|&x| *x == true).count();
        let m_off = terms.iter().filter(|&x| *x == false).count();
        Terms {terms, m, m_on, m_off}

    }
}
impl State for Terms
{
    type Move = ToggleSwap;

    fn draw_move<R : Rng>(&self, rng : &mut R) -> ToggleSwap {
        // Every possible state transition is equally likely.
        let n = self.m + self.m_on * self.m_off;
        let x = rng.random_range(0..n);

        // In m cases we flip the state of a single term.
        if x < self.m {
            Toggle(x)
        }
        // in m_on * m_off cases we swap the states of two terms with different states.
        else {
            // map random number x to pairs of indices a, b
            let k = x - self.m;

            let (i, j) = find_swap_indices(k, self.m_on, self.m_off, &self.terms)
                .expect("Swap index out of range");

            Swap(i, j)
        }
    }

    fn apply(&mut self, m: & ToggleSwap) {
        match *m {
            Toggle(i) => {
                self.terms[i] = !self.terms[i]
            }
            Swap(i_on, i_off) => {
                self.terms[i_on] = false;
                self.terms[i_off] = true;
            }
        }
    }

    fn revert(&mut self, m: & ToggleSwap) {
        match *m {
            Toggle(i) => {
                self.terms[i] = !self.terms[i]
            }
            Swap(i_on, i_off) => {
                self.terms[i_on] = true;
                self.terms[i_off] = false;
            }
        }
    }
}


/// Maps a random number `k` to indices (index_on, index_off) in the terms vector.
/// Assumption: k < m_on * m_off
fn find_swap_indices(k : usize, m_on : usize, m_off : usize, terms: &[bool]) -> Option<(usize, usize)> {
    if k >= m_on * m_off {
        return None;
    }

    let target_on_nth = k % m_on;
    let target_off_nth = k / m_on;

    let mut on_idx = None;
    let mut off_idx = None;
    let mut current_on_count = 0;
    let mut current_off_count = 0;

    // Scan the vector once to find the actual indices of "nth active" and "nth inactive" terms.
    for (i, &is_active) in terms.iter().enumerate() {
        if is_active {
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

    #[test]
    fn test_swap_indices_mapping() {
        let terms = vec![true, false, false, true, false];
        let m_on = 2;
        let m_off = 3;

        // Case k=0: target_on=0, target_off=0 -> Expect indices (0, 1)
        assert_eq!(find_swap_indices(0, m_on, m_off, &terms).unwrap(), (0, 1));

        // Case k=1: target_on=1, target_off=0 -> Expect indices (3, 1)
        assert_eq!(find_swap_indices(1, m_on, m_off, &terms).unwrap(), (3, 1));

        // Case k=2: target_on=0, target_off=1 -> Expect indices (0, 2)
        assert_eq!(find_swap_indices(2, m_on, m_off, &terms).unwrap(), (0, 2));

        // Case k=3: target_on=1, target_off=1 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(3, m_on, m_off, &terms).unwrap(), (3, 2));

        // Case k=4: target_on=0, target_off=2 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(4, m_on, m_off, &terms).unwrap(), (0, 4));

        // Case k=5: target_on=1, target_off=2 -> Expect indices (3, 2)
        assert_eq!(find_swap_indices(5, m_on, m_off, &terms).unwrap(), (3, 4));
    }

    #[test]
    #[should_panic]
    fn test_out_of_bounds() {
        // If we ask for a k that implies more active terms than exist
        let terms = vec![true, false];
        let m_on = 1;
        let m_off = 1;

        find_swap_indices(5, m_on, m_off, &terms).unwrap();
    }
}