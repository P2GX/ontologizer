use oboannotation::go::GoAnnotations;
use oboannotation::go::stats::get_annotation_map;
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};

pub struct GeneSet {
    recognized_genes: HashSet<String>,
    unrecognized_genes: HashSet<String>,
}

impl GeneSet {
    pub fn from_file(path: &str, annotations: &GoAnnotations) -> Result<Self, String> {
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

        // Separate Genes by GoAnnotations
        let annotations = get_annotation_map(&annotations);

        let (recognized_genes, unrecognized_genes): (HashSet<String>, HashSet<String>) =
            genes.into_iter().partition(|g| annotations.contains_key(g));

        Ok(GeneSet {
            recognized_genes,
            unrecognized_genes,
        })
    }

    pub fn recognized_genes(&self) -> &HashSet<String> {
        &self.recognized_genes
    }

    pub fn unrecognized_genes(&self) -> &HashSet<String> {
        &self.unrecognized_genes
    }
}
