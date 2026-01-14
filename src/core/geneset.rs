use oboannotation::{
    go::{GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

// Tuple struct for gene symbols that ensure normalization (uppercase, trimmed, no comments)
#[derive(Debug)]
pub struct GeneSet {
    recognized: HashSet<String>,
    unrecognized: HashSet<String>,
}

impl GeneSet {
    pub fn recognized_genes(&self) -> &HashSet<String> {
        &self.recognized
    }
    pub fn unrecognized_genes(&self) -> &HashSet<String> {
        &self.unrecognized
    }
    pub fn gene_count(&self) -> usize {
        self.recognized.len() + self.unrecognized.len()
    }
}

// Loads gene symbols from a text file. Each line in the file should contain a single gene symbol.
pub fn load_gene_set(path: &str) -> Result<HashSet<String>, String> {
    let file =
        File::open(path).map_err(|err| format!("Failed to open file '{}': {}", path, err))?;
    let reader = BufReader::new(file);

    let mut genes: HashSet<String> = HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|err| format!("Failed to read line {}: {}", i + 1, err))?;
        let gene = line.trim().to_string();
        if !gene.is_empty() {
            genes.insert(gene);
        }
    }
    Ok(genes)
}

pub fn separate_gene_set(
    annotated_genes: &HashSet<String>,
    genes: HashSet<String>,
) -> GeneSet {
    let mut recognized: HashSet<String> = HashSet::new();
    let mut unrecognized: HashSet<String> = HashSet::new();

    for sym in genes {
        if annotated_genes.contains(sym.as_str()) {
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

        let annotated_genes = get_annotation_map(&annotations).into_keys().collect();

        let genes = load_gene_set(gene_set_path).expect("Failed to load gene set");
        let gene_set = separate_gene_set(&annotated_genes, genes);

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
