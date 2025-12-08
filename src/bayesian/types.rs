use crate::core::{AnnotationIndex, GeneSymbol};
use anyhow::Result;
use ontolius::TermId;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

/// Fixed MGSA model parameters.
pub struct MgsaParams {
    pub alpha: f64,
    pub beta: f64,
    pub q: f64,
}

/// Fixed problem definition. Only place that knows about AnnotationIndex.
pub struct Problem<'a> {
    pub ann: &'a AnnotationIndex,
    pub genes: HashSet<GeneSymbol>,
}

impl<'a> Problem<'a> {
    pub fn all_genes(&self) -> HashSet<GeneSymbol> {
        self.ann.genes_to_terms.keys().cloned().collect()
    }

    pub fn all_terms(&self) -> HashSet<TermId> {
        self.ann.terms_to_genes.keys().cloned().collect()
    }

    pub fn observed_genes(&self) -> HashSet<GeneSymbol> {
        self.genes.clone()
    }
}

/// MCMC configuration parameters.
pub struct MgsaConfig {
    pub steps: u64,
    pub burn_in: u64,
    // todo!(should allow for a seed)
}

/// Variable state during sampling: term on/off.
/// Hidden gene states are derived on the fly when needed.
#[derive(Clone)]
pub struct State {
    pub terms: HashMap<TermId, bool>, // key set fixed after init; values flip
    pub hidden: HashMap<GeneSymbol, bool>,
}

impl State {
    pub fn new(terms: &HashSet<TermId>) -> Self {
        let terms = terms.iter().cloned().map(|t| (t, false)).collect();
        let hidden = HashMap::new();
        Self { terms, hidden }
    }
}

/// Accumulator of term counts across kept samples.
pub struct MgsaResult {
    pub counts: HashMap<TermId, usize>,
    pub samples: usize,
}

impl MgsaResult {
    pub fn new(terms: &HashSet<TermId>) -> Self {
        let counts = terms.iter().cloned().map(|t| (t, 0usize)).collect();
        Self { counts, samples: 0 }
    }

    pub fn update(&mut self, state: &State) {
        for (term, &on) in state.terms.iter() {
            if on {
                // safe: keys identical to initialization
                *self.counts.get_mut(term).unwrap() += 1;
            }
        }
        self.samples += 1;
    }

    pub fn write_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let start = Instant::now();

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Collect and sort by posterior desc
        let mut rows: Vec<(&TermId, usize, f64)> = self
            .counts
            .iter()
            .map(|(t, c)| {
                let post = if self.samples > 0 {
                    *c as f64 / self.samples as f64
                } else {
                    f64::NAN
                };
                (t, *c, post)
            })
            .collect();

        rows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

        // Header
        writeln!(writer, "term_id\tcount\tsamples\tposterior")?;

        // Rows
        for (term, count, post) in rows {
            if post.is_finite() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{:.6e}",
                    term, count, self.samples, post
                )?;
            } else {
                // samples == 0 case
                writeln!(writer, "{}\t{}\t{}\tNA", term, count, self.samples)?;
            }
        }

        eprintln!("MGSA: wrote results in {:?}", start.elapsed());
        Ok(())
    }
}
