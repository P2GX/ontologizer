mod annotations;
mod util;
pub(crate) mod result;
mod problem;

use std::collections::HashSet;
pub use util::{load_gene_set, overlap_sets};

pub use annotations::AnnotationIndex;