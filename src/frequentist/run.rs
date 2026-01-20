use indexmap::IndexSet;
use std::collections::HashSet;

use crate::core::AnnotationIndex;

use crate::core::result::EnrichmentResult;
use crate::frequentist::algorithm::{OneSidedEnrichmentTest, StatisticalTest};
use crate::frequentist::correction::{BonferroniHolm, Correction};
use crate::frequentist::distribution::{Hypergeometric, LogFactorialCache};
use crate::frequentist::measure::PValue;
use ontolius::ontology::csr::FullCsrOntology;

#[allow(non_snake_case)]
pub fn run(
    ontology: &FullCsrOntology,
    annotations: AnnotationIndex,
    study_genes: HashSet<String>,
) -> EnrichmentResult {
    let n = study_genes.len();
    let N = annotations.get_genes().len();

    let cache = LogFactorialCache::new(N);
    let test = OneSidedEnrichmentTest;

    let observed_genes: Vec<bool> = (0..N)
        .map(|i| study_genes.contains(annotations.get_gene_by_index(i)))
        .collect();

    let study_gene_indices: IndexSet<usize> = study_genes
        .iter()
        .filter_map(|gene| annotations.get_index_by_gene(gene))
        .collect();

    let terms_to_genes = annotations.get_terms_to_genes(true);

    let mut measures: Vec<PValue> = Vec::new();

    for population_genes_indices in terms_to_genes {
        let K = population_genes_indices.len();

        let k = population_genes_indices
            .intersection(&study_gene_indices)
            .count();

        let hypergeo = Hypergeometric::new(n, K, N, &cache);
        let raw_p_value = test.calculate(&hypergeo, k);

        let counts = (k as u32, n as u32, K as u32, N as u32);
        let measure = PValue::new(raw_p_value, counts);

        measures.push(measure);
    }

    let correction = BonferroniHolm;
    correction.adjust(&mut measures, |m| &mut m.pvalue);

    // Create the Result
    let mut result =
        EnrichmentResult::from_measure(&measures, &ontology, &annotations, &observed_genes);

    // Optional: Sort
    result.sort_by_score(false); // ascending for p-value

    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::core;
    use csv::Writer;
    use oboannotation::go::stats::get_annotation_map;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use std::process;

    #[test]
    fn test_go_analysis() {
        // Paths to the necessary files
        // These paths should be adjusted based on your local setup
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let population_genes_path = "tests/data/GOnone/population.txt";
        let study_genes_path = "tests/data/GOnone/study.txt";

        // ------ Load Gene Sets ------
        let study_genes = core::load_gene_set(study_genes_path).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        let population_genes = core::load_gene_set(population_genes_path).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

        // ------ Load Gene Ontology and Annotations ------
        let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
        let ontology: FullCsrOntology =
            ontology_loader
                .load_from_path(go_path)
                .unwrap_or_else(|err| {
                    eprintln!("Error: {}", err);
                    process::exit(1);
                });

        let annotations_loader = GoGafAnnotationLoader;
        let annotations: GoAnnotations = annotations_loader
            .load_from_path(gaf_path)
            .expect("Could not load GAF file");

        let annotation_index =
            AnnotationIndex::new(annotations, &ontology, Some(&population_genes));
        let annotated_genes = get_annotation_map(&annotation_index.annotations)
            .into_keys()
            .collect();

        let study_gene_set = core::overlap_sets(&annotated_genes, &study_genes);

        let result = run(&ontology, annotation_index, study_gene_set);

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GOnone/results_f.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }
}
