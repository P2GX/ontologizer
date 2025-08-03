use std::collections::HashSet;
use std::time::Instant;

use crate::{annotations::AnnotationContainer, gene_set::GeneSet, hypergeom::Hypergeometric};

use ontolius::{
    ontology::{OntologyTerms, csr::FullCsrOntology},
    term::MinimalTerm,
};

use super::pvalue_calculation::PValueCalculation;
use super::results::{AnalysisResults, GOTermResult, MethodEnum, MtcEnum};

pub struct TermForTerm;
impl PValueCalculation for TermForTerm {
    fn calculate_p_values(
        &self,
        go: &FullCsrOntology,
        annotation_container: &AnnotationContainer,
        study: &GeneSet,
        population: &GeneSet,
    ) -> AnalysisResults {
        let start = Instant::now();
        let mut results = AnalysisResults::new(MethodEnum::TermForTerm, MtcEnum::None); // Placeholder for MTC method
        let study_genes = study.gene_symbols();
        let pop_genes = population.gene_symbols();
        let n = study_genes.len();
        let m = pop_genes.len();

        for (term, nt) in annotation_container.study_annotations().iter() {
            if *nt > 1 {
                let annotated_genes = annotation_container.term_genes_map().get(term).unwrap();

                let nt = gene_set_intersection(study_genes, &annotated_genes).len(); // Eigentlich ist die Funktion (Schnittmenge) überflüssig, da die study_annotations map schon die nt Werte enthält
                let mt = annotated_genes.len();
                let mut hypergeom = Hypergeometric::new();
                let raw_p_value = hypergeom.phyper(nt - 1, m, mt, n, false).unwrap();

                let term_id = go.term_by_id(term);
                let result = GOTermResult::new(
                    term_id.unwrap().name().to_string(),
                    term.to_string(),
                    (nt as u32, mt as u32, n as u32, m as u32),
                    raw_p_value as f32,
                    raw_p_value as f32, // use raw pvalue as default
                );

                results.add_result(result);
            }
        }
        results.sort_by_p_value();
        let duration = start.elapsed();
        eprintln!("P-value calculation took: {:?}", duration);

        results
    }
}
// Computes the intersection of a gene set with the set of genes annotated to a specific GO term.
fn gene_set_intersection(
    gene_set: &HashSet<String>,
    genes_annotated_to_term: &HashSet<String>,
) -> HashSet<String> {
    let intersection: HashSet<_> = gene_set
        .intersection(genes_annotated_to_term)
        .cloned()
        .collect();
    intersection
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::{
        gene_set::{build_gene_set, load_gene_file},
        ontologizer::Ontologizer,
    };

    #[test]
    #[ignore]
    fn test_go_analysis() {
        // Paths to the necessary files
        // These paths should be adjusted based on your local setup
        let go_path = "tests/data/go-basic.json";
        let gaf_path = "tests/data/goa_human.gaf";
        let pop_set_path = "tests/data/population.txt";
        let study_set_path = "tests/data/study.txt";

        // Load the GO ontology
        let go = Ontologizer::new(go_path);
        let go_ref = go.ontology();

        // Load the GOA annotations
        let mut annotation_container = AnnotationContainer::new(gaf_path);

        // Load the population and study gene sets
        let study_set_list = load_gene_file(study_set_path).expect("Failed to load study list");
        let study_set = build_gene_set(study_set_list, &annotation_container.annotations());
        let pop_set_list = load_gene_file(pop_set_path).expect("Failed to load population list");
        let pop_set = build_gene_set(pop_set_list, &annotation_container.annotations());

        // Build map that contains all GO terms annotated in the study set and their counts.
        // (we only want to analyze terms that are annotated in the study set)
        annotation_container.build_study_annotations(&study_set, go_ref);
        // Build map that contains all GO terms annotated in the population set and their associated genes (in population set).
        annotation_container.build_term_genes_map(&pop_set, go_ref);

        

        let analysis_method = TermForTerm;
        let results =
            analysis_method.calculate_p_values(go_ref, &annotation_container, &study_set, &pop_set);
        results
            .write_tsv("tests/data/enrichment_results.tsv")
            .expect("Failed to write results");
    }
}
