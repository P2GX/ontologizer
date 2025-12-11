pub trait Observation{

}

pub(crate) struct Genes{
    // Computational vector used by *ALGORITHM*. Mapping gene indices to names by *ANNOTATIONS**
    genes : Vec<bool>,
    n: usize,
    n_on: usize,
    n_off: usize,
}

impl Observation for Genes{
    
}