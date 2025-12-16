use crate::core::{AnnotationIndex, GeneSymbol, Ontologizer, load_gene_set, separate_gene_set};
use oboannotation::go::GoGafAnnotationLoader;
use oboannotation::go::stats::get_annotation_map;
use oboannotation::io::AnnotationLoader;

struct Result {
    term_probs: Vec<f32>,
    term_genes: Vec<GeneSymbol>,
}

struct MgsaParameter {
    p: f32,
    alpha: f32,
    beta: f32,
}

pub fn run_mgsa(parameter: MgsaParameter) -> Result {
    let result: Result;

    let go_path = "tests/data/go-basic.json";
    let gaf_path = "tests/data/goa_human.gaf";
    let pop_set_path = "tests/data/population.txt";
    let study_set_path = "tests/data/study.txt";

    // Load the GO ontology
    let go = Ontologizer::new(go_path);
    let go_ref = go.ontology();

    // Load the GOA annotations
    let goa_path = "tests/data/goa_human.gaf";
    let annotations = GoGafAnnotationLoader
        .load_from_path(goa_path)
        .expect("Could not load GAF file");

    let annotated_genes = get_annotation_map(&annotations).into_keys().collect();

    // Load the population and study gene sets
    let gene_symbols = load_gene_set(study_set_path).expect("Failed to parse study gene set");
    let genes = separate_gene_set(&annotated_genes, gene_symbols);

    todo!()
}
