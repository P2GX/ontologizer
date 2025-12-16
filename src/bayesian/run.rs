use crate::core::{load_gene_set, separate_gene_set, AnnotationIndex, GeneSymbol, Ontologizer};

struct Result{
    term_probs : Vec<f32>,
    term_genes : Vec<GeneSymbol>
}

struct MgsaParameter{
    p : f32,
    alpha : f32,
    beta : f32
}

pub fn run_mgsa(parameter : MgsaParameter) -> Result {
    let result : Result;

    let go_path = "tests/data/go-basic.json";
    let gaf_path = "tests/data/goa_human.gaf";
    let pop_set_path = "tests/data/population.txt";
    let study_set_path = "tests/data/study.txt";

    // Load the GO ontology
    let go = Ontologizer::new(go_path);
    let go_ref = go.ontology();

    // Load the GOA annotations
    let mut annotations = AnnotationIndex::new(gaf_path, go_ref);

    // Load the population and study gene sets
    let gene_symbols =
        load_gene_set(study_set_path).expect("Failed to parse study gene set");
    let genes =
        separate_gene_set(&annotations.annotations, gene_symbols);


    let terms = annotations.terms();

    todo!()
}