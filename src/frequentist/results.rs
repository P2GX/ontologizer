use std::{
    fmt,
    fs::File,
    io::{BufWriter, Result, Write},
    time::Instant,
};

use ontolius::{
    TermId,
    common::go::{BIOLOGICAL_PROCESS, CELLULAR_COMPONENT, MOLECULAR_FUNCTION},
    ontology::{HierarchyQueries, csr::FullCsrOntology},
};
use serde::Serialize;

// This module provides the results of enrichment analysis, including methods to write results to a TSV file.
// Every entry in the results corresponds to one GO term
#[derive(Serialize, Clone)]
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

    pub fn results(&self) -> &Vec<GOTermResult> {
        &self.results
    }

    pub fn mtc_method(&self) -> &MtcEnum {
        &self.mtc_method
    }

    pub fn analysis_method(&self) -> &MethodEnum {
        &self.analysis_method
    }

    pub fn num_significant_results(&self) -> usize {
        self.results.iter().filter(|r| r.adj_pval <= 0.05).count()
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
                // result.p_val(),
                // result.adj_pval(),
            )?;
        }
        let duration = start.elapsed();
        eprintln!("Writing results to {} took: {:?}", path, duration);

        Ok(())
    }
}

// Represents a single analysis result, including the Analysis and Multiple Testing Correction method used,
// term information, counts, and p-values.
#[derive(Serialize, Clone)]
pub struct GOTermResult {
    term_label: String,
    term_id: String,
    counts: (u32, u32, u32, u32), // nt, mt, n, m (annotated genes in study, annotated genes in population, total genes in study, total genes in population)
    p_val: f32,                   // raw p-value
    adj_pval: f32,
    aspect: &'static str, // Aspect of the GO term: BP, MF, CC
}

impl GOTermResult {
    pub fn new(
        term_label: String,
        term_id: String,
        counts: (u32, u32, u32, u32),
        p_val: f32,
        adj_pval: f32,
        aspect: &'static str,
    ) -> Self {
        Self {
            term_label,
            term_id,
            counts,
            p_val,
            adj_pval,
            aspect,
        }
    }

    pub fn set_adj_pval(&mut self, val: f32) {
        self.adj_pval = val;
    }

    pub fn p_val(&self) -> &f32 {
        &self.p_val
    }

    pub fn adj_pval(&self) -> f32 {
        self.adj_pval
    }

    pub fn term_id(&self) -> &str {
        &self.term_id
    }

    pub fn term_label(&self) -> &str {
        &self.term_label
    }
}

pub fn get_term_aspect(go: &FullCsrOntology, term: &TermId) -> &'static str {
    if go.is_equal_or_descendant_of(term, &BIOLOGICAL_PROCESS) {
        "BP"
    } else if go.is_equal_or_descendant_of(term, &MOLECULAR_FUNCTION) {
        "MF"
    } else if go.is_equal_or_descendant_of(term, &CELLULAR_COMPONENT) {
        "CC"
    } else {
        "Unknown"
    }
}

// Represents the method used for the analysis, such as Term for Term or Parent Child.
#[derive(Serialize, Debug, Clone)]
pub enum MethodEnum {
    TermForTerm,
    ParentChildUnion,
    ParentChildIntersection,
    MGSA,
}

// Represents the method used for multiple testing correction, such as Bonferroni or FDR.
#[derive(Serialize, Debug, Clone)]
pub enum MtcEnum {
    Bonferroni,
    BonferroniHolm,
    BenjaminiHochberg,
    None, // No multiple testing correction
}

impl MethodEnum {
    pub fn new(name: String) -> Self {
        match name.as_str() {
            "TermForTerm" => MethodEnum::TermForTerm,
            "ParentChildUnion" => MethodEnum::ParentChildUnion,
            "ParentChildIntersection" => MethodEnum::ParentChildIntersection,
            "MGSA" => MethodEnum::MGSA,
            _ => panic!("Unknown MethodEnum variant: {}", name),
        }
    }
}

impl MtcEnum {
    pub fn new(name: String) -> Self {
        match name.as_str() {
            "Bonferroni" => MtcEnum::Bonferroni,
            "BonferroniHolm" => MtcEnum::BonferroniHolm,
            "BenjaminiHochberg" => MtcEnum::BenjaminiHochberg,
            "None" => MtcEnum::None,
            _ => panic!("Unknown MtcEnum variant: {}", name),
        }
    }
}

impl fmt::Display for MethodEnum {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MethodEnum::TermForTerm => "TermForTerm",
            MethodEnum::ParentChildUnion => "ParentChildUnion",
            MethodEnum::ParentChildIntersection => "ParentChildIntersection",
            MethodEnum::MGSA => "MGSA",
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
            MtcEnum::None => "None",
        };
        write!(f, "{}", s)
    }
}
