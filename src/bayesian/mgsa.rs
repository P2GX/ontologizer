

#[cfg(test)]
mod test {
    use crate::core::AnnotationIndex;
    use crate::core::{load_gene_set, separate_gene_set};
    use crate::core::Ontologizer;
    use super::*;


    #[test]
    fn test_mgsa() {
        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let pop_set_path = "tests/data/population.txt";
        let study_set_path = "tests/data/study.txt";

        let gene_ontology = Ontologizer::new(go_path);

        // Load the GOA annotations
        let mut annotations = AnnotationIndex::new(gaf_path);

        // Load the population and study gene sets

        // Load the population and study gene sets
        let study_gene_symbols = load_gene_set(study_set_path).expect("Failed to parse study gene set");
        let study_gene_set = separate_gene_set(&annotations.annotations(), study_gene_symbols);

        let pop_gene_symbols = load_gene_set(pop_set_path).expect("Failed to parse population gene set");
        let pop_gene_set = separate_gene_set(&annotations.annotations(), pop_gene_symbols);

        // Build map that contains all GO terms annotated in the study set and their counts.
        // (we only want to analyze terms that are annotated in the study set)
        annotations.compute_term_counts(&study_gene_set, gene_ontology.ontology());
        // Build map that contains all GO terms annotated in the population set and their associated genes (in population set).
        annotations.build_terms_to_genes(&pop_gene_set, gene_ontology.ontology());

        // Mgsa.calculate_probabilities(gene_ontology.ontology(), &annotations, &study_gene_set, &pop_gene_set)
    }
}