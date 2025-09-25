use oboannotation::{
    go::{GoAnnotations, GoGafAnnotationLoader, stats::get_annotation_map},
    io::{AnnotationLoadError, AnnotationLoader},
};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

#[derive(Debug)]
/// A struct representing a set of genes (population set, study set), including recognized and unrecognized gene symbols.
pub struct GeneSet {
    gene_symbols: HashSet<String>,
    unrecognized_gene_symbols: HashSet<String>,
    del_duplicates: usize,
}

impl GeneSet {
    pub fn gene_symbols(&self) -> &HashSet<String> {
        &self.gene_symbols
    }
    pub fn unrecognized_gene_symbols(&self) -> &HashSet<String> {
        &self.unrecognized_gene_symbols
    }

    pub fn gene_count(&self) -> usize {
        self.gene_symbols.len() + self.unrecognized_gene_symbols.len()
    }

    pub fn del_duplicates(&self) -> usize {
        self.del_duplicates
    }
}

// Loads gene symbols from a text file.
// Each line in the file should contain a single gene symbol.
pub fn load_gene_file(path: &str) -> Result<Vec<String>, String> {
    let file =
        File::open(path).map_err(|err| format!("Failed to open file '{}': {}", path, err))?;
    let reader = BufReader::new(file);
    let mut genes = vec![];

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|err| format!("Failed to read line {}: {}", i + 1, err))?;
        let gene = line.trim().to_string();
        if !gene.is_empty() {
            genes.push(gene);
        }
    }
    Ok(genes)
}

// Builds a 'GeneSet' by separating recognized and unrecognized gene symbols
// based on the provided GO annotations.
// Genes that are present in the annotation map are added to 'gene_symbols',
// while genes that are not found in the annotation map are added to 'unrecognized_gene_symbols'.
pub fn build_gene_set(gene_list: Vec<String>, annotations: &GoAnnotations) -> GeneSet {
    let annotation_map = get_annotation_map(&annotations);

    let mut gene_symbols = HashSet::new();
    let mut unrecognized_gene_symbols = HashSet::new();
    let mut deleted_duplicates: usize = 0;

    for gene in gene_list {
        if annotation_map.contains_key(&gene) {
            if !gene_symbols.insert(gene) {
                deleted_duplicates += 1;
            }
        } else {
             if !unrecognized_gene_symbols.insert(gene) {
                deleted_duplicates += 1;
            }
        }
    }
    GeneSet {
        gene_symbols: gene_symbols,
        unrecognized_gene_symbols: unrecognized_gene_symbols,
        del_duplicates: deleted_duplicates,
    }
}

// Loads GOA annotations from a .gaf file and creates a `GoAnnotations` object.
// Currently only supports uncompressed GOA files.
pub fn load_goa_annotations(gaf_path: &str) -> Result<GoAnnotations, AnnotationLoadError> {
    let loader = GoGafAnnotationLoader;
    let annotations = loader.load_from_path(gaf_path)?;
    Ok(annotations)
}

#[cfg(test)]
mod test {
    use super::*;
    use rstest::rstest;

    #[test]
    fn test_load_gene_file() {
        let gene_list_path = "tests/data/study.txt";
        let genes = load_gene_file(gene_list_path).unwrap();
        eprintln!("Loaded {} genes from {}", genes.len(), gene_list_path);
        assert!(!genes.is_empty(), "Gene list should not be empty");
    }

    #[test]
    fn test_load_goa_annotations() {
        let gaf_path = "tests/data/goa_human.gaf";

        let annotations = load_goa_annotations(gaf_path).expect("GAF could not be loaded");

        eprintln!(
            "Loaded {} GOA annotations from {}",
            annotations.annotations.len(),
            gaf_path
        );

        assert_eq!(
            !annotations.annotations.is_empty(),
            true,
            "GOA annotations should not be empty"
        );

        // eprintln!("First annotation: {:#?}", annotations.annotations.first());
    }

    #[rstest]
    #[case::studyfile("tests/data/study.txt")]
    #[case::populationfile("tests/data/population.txt")] // file only contains Ensembl IDs -> unrecognized genes
    fn test_build_gene_set(#[case] gene_set_path: &str) {
        let goa_path = "tests/data/goa_human.gaf";
        let annotations = load_goa_annotations(goa_path).expect("Failed to load GOA annotations");

        let gene_list = load_gene_file(gene_set_path).expect("Failed to load study gene list");
        let gene_set = build_gene_set(gene_list, &annotations);

        eprintln!(
            "Built gene set with {} recognized genes, {} unrecognized genes and {} deleted duplicates",
            gene_set.gene_symbols.len(),
            gene_set.unrecognized_gene_symbols.len(),
            gene_set.del_duplicates()

        );

        assert!(
            !gene_set.gene_symbols.is_empty() || !gene_set.unrecognized_gene_symbols.is_empty(),
            "Gene set should not be empty in {:?}",
            gene_set_path
        );
    }
}
