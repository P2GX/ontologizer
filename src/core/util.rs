use oboannotation::{
    go::{GoGafAnnotationLoader, stats::get_annotation_map},
    io::AnnotationLoader,
};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
};
use std::hash::Hash;

pub fn overlap_sets<T>(set_a : &HashSet<T>, set_b: &HashSet<T>) -> HashSet<T>
where T: Eq + Hash + Clone
{
    set_a.intersection(set_b).cloned().collect()
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
