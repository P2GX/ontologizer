use oboannotation::{
    go::{GoAnnotations, GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

// Tuple struct for gene symbols that ensure normalization (uppercase, trimmed, no comments)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]

pub struct GeneSymbol(String);

impl GeneSymbol {
    // Ensures that gene symbols do not start with '#' and have no leading or trailing whitespaces.
    pub fn new<S: AsRef<str>>(s: S) -> Option<Self> {
        let s = s.as_ref().trim_start();
        if s.is_empty() || s.starts_with('#') {
            return None;
        }
        let s = s.trim();
        let mut out = s.to_owned();
        Some(Self(out))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for GeneSymbol {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

#[derive(Debug)]
pub struct GeneSet {
    recognized: HashSet<GeneSymbol>,
    unrecognized: HashSet<GeneSymbol>,
}

impl GeneSet {
    pub fn recognized_genes(&self) -> &HashSet<GeneSymbol> {
        &self.recognized
    }
    pub fn unrecognized_genes(&self) -> &HashSet<GeneSymbol> {
        &self.unrecognized
    }
    pub fn gene_count(&self) -> usize {
        self.recognized.len() + self.unrecognized.len()
    }
}

// Loads gene symbols from a text file. Each line in the file should contain a single gene symbol.
pub fn load_gene_set(path: &str) -> Result<HashSet<GeneSymbol>, String> {
    let file =
        File::open(path).map_err(|err| format!("Failed to open file '{}': {}", path, err))?;
    let reader = BufReader::new(file);

    let mut genes: HashSet<GeneSymbol> = HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|err| format!("Failed to read line {}: {}", i + 1, err))?;

        // Normalize/filter (skips blanks/comments)
        match GeneSymbol::new(&line) {
            Some(sym) => genes.insert(sym),
            None => continue,
        };
    }
    Ok(genes)
}

pub fn separate_gene_set(annotations: &GoAnnotations, genes: HashSet<GeneSymbol>) -> GeneSet {
    let mut recognized: HashSet<GeneSymbol> = HashSet::new();
    let mut unrecognized: HashSet<GeneSymbol> = HashSet::new();

    let annotation_map = get_annotation_map(annotations);

    for sym in genes {
        if annotation_map.contains_key(sym.as_str()) {
            // Insert; if already present, it's a duplicate
            recognized.insert(sym);
        } else {
            unrecognized.insert(sym);
        }
    }
    GeneSet {
        recognized,
        unrecognized,
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::studyfile("tests/data/study.txt")]
    #[case::populationfile("tests/data/population.txt")] // file only contains Ensembl IDs -> unrecognized genes
    fn test_build_gene_set(#[case] gene_set_path: &str) {
        let goa_path = "tests/data/goa_human.gaf";
        let annotations = GoGafAnnotationLoader
            .load_from_path(goa_path)
            .expect("Could not load GAF file");

        let genes = load_gene_set(gene_set_path).expect("Failed to load gene set");
        let gene_set = separate_gene_set(&annotations, genes);

        eprintln!(
            "Built gene set with {} recognized genes, {} unrecognized genes.",
            gene_set.recognized.len(),
            gene_set.unrecognized.len(),
        );

        assert!(
            !gene_set.recognized.is_empty() || !gene_set.unrecognized.is_empty(),
            "Gene set should not be empty in {:?}",
            gene_set_path
        );
    }
}
