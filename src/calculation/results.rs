use std::{
    fmt,
    fs::File,
    io::{BufWriter, Result, Write},
    time::Instant,
};

// This module provides the results of enrichment analysis, including methods to write results to a TSV file.
// Every entry in the results corresponds to one GO term
pub struct AnalysisResults {
    results: Vec<GOTermResult>,
    mtc_method: MtcEnum,         // Method used for multiple testing correction
    analysis_method: MethodEnum, // Method used for the analysis
}

impl AnalysisResults {
    pub fn new(analysis_method: MethodEnum, mtc_method: MtcEnum) -> Self {
        Self {
            results: Vec::new(),
            mtc_method,
            analysis_method,
        }
    }

    pub fn add_result(&mut self, result: GOTermResult) {
        self.results.push(result);
    }

    pub fn num_hypotheses(&self) -> f32 {
        self.results.len() as f32
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<GOTermResult> {
        self.results.iter_mut()
    }

    pub fn sort_by_p_value(&mut self) {
        self.results
            .sort_by(|a, b| a.p_val.partial_cmp(&b.p_val).unwrap());
    }
    // Writes the enrichment results to a TSV file.
    pub fn write_tsv(&self, path: &str) -> Result<()> {
        let start = Instant::now();
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Header
        writeln!(
            writer,
            "analysis_method\tmtc_method\tterm_label\tterm_id\tnt\tmt\tn\tm\tp_val\tadj_pval"
        )?;

        for result in &self.results {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.analysis_method,
                self.mtc_method,
                result.term_label,
                result.term_id,
                result.counts.0,
                result.counts.1,
                result.counts.2,
                result.counts.3,
                format!("{:e}", result.p_val),
                format!("{:e}", result.adj_pval),
            )?;
        }
        let duration = start.elapsed();
        eprintln!("Writing results to {} took: {:?}", path, duration);

        Ok(())
    }
}

// Represents a single analysis result, including the Analysis and Multiple Testing Correction method used,
// term information, counts, and p-values.
pub struct GOTermResult {
    term_label: String,
    term_id: String,
    counts: (u32, u32, u32, u32), // nt, mt, n, m (annotated genes in study, annotated genes in population, total genes in study, total genes in population)
    p_val: f32,                   // raw p-value
    adj_pval: f32,
}

impl GOTermResult {
    pub fn new(
        term_label: String,
        term_id: String,
        counts: (u32, u32, u32, u32),
        p_val: f32,
        adj_pval: f32,
    ) -> Self {
        Self {
            term_label,
            term_id,
            counts,
            p_val,
            adj_pval,
        }
    }

    pub fn set_adj_pval(&mut self, val: f32) {
        self.adj_pval = val;
    }

    pub fn p_val(&self) -> &f32 {
        &self.p_val
    }
}

// Represents the method used for the analysis, such as Term for Term or Parent Child.
pub enum MethodEnum {
    TermForTerm,
    ParentChild,
}

// Represents the method used for multiple testing correction, such as Bonferroni or FDR.
pub enum MtcEnum {
    Bonferroni,
    BonferroniHolm,
    BenjaminiHochberg,
    FDR,
    None, // No multiple testing correction
}

impl fmt::Display for MethodEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MethodEnum::TermForTerm => "TermForTerm",
            MethodEnum::ParentChild => "ParentChild",
        };
        write!(f, "{}", s)
    }
}

impl fmt::Display for MtcEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MtcEnum::Bonferroni => "Bonferroni",
            MtcEnum::BonferroniHolm => "BonferroniHolm",
            MtcEnum::BenjaminiHochberg => "BenjaminiHochberg",
            MtcEnum::FDR => "FDR",
            MtcEnum::None => "None",
        };
        write!(f, "{}", s)
    }
}
