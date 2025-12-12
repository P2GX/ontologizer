pub trait Observation {
    type Value: Copy;
    fn get(&self, i: usize) -> Self::Value;
}

pub(crate) struct Genes {
    // Computational vector used by *ALGORITHM*. Mapping gene indices to names by *ANNOTATIONS**
    genes: Vec<bool>,
    n: usize,
    n_on: usize,
    n_off: usize,
}

impl Genes {
    pub(crate) fn new(genes: Vec<bool>) -> Genes {
        let n = genes.len();
        let n_on = genes.iter().filter(|&x| *x == true).count();
        let n_off = genes.iter().filter(|&x| *x == false).count();
        Genes {
            genes,
            n,
            n_on,
            n_off,
        }
    }
}

impl Observation for Genes {
    type Value = bool;

    fn get(&self, i: usize) -> bool {
        self.genes[i]
    }
}
