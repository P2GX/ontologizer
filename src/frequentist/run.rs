use indexmap::IndexSet;
use std::collections::HashSet;

use crate::core::AnnotationIndex;

use crate::core::result::EnrichmentResult;
use crate::frequentist::algorithm::{OneSidedEnrichmentTest, StatisticalTest};
use crate::frequentist::correction;
use crate::frequentist::correction::Correction;
use crate::frequentist::distribution::{Hypergeometric, LogFactorialCache};
use crate::frequentist::measure::PValue;
use ontolius::ontology::csr::FullCsrOntology;

#[allow(non_snake_case)]
pub fn analysis(
    ontology: &FullCsrOntology,
    annotations: AnnotationIndex,
    study_genes: HashSet<String>,
) -> EnrichmentResult {
    let n = study_genes.len();
    let N = annotations.get_genes().len();

    let cache = LogFactorialCache::new(N);
    let test = OneSidedEnrichmentTest;

    let observed_genes: Vec<bool> = (0..N)
        .map(|i| study_genes.contains(annotations.get_gene_by_idx(i)))
        .collect();

    let study_gene_indices: IndexSet<usize> = study_genes
        .iter()
        .filter_map(|gene| annotations.get_idx_by_gene(gene))
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

    let correction = correction::None;
    correction.adjust(&mut measures, |m| &mut m.pvalue);

    // Create the Result
    let mut result =
        EnrichmentResult::from_measures(&measures, &ontology, &annotations, &observed_genes);

    // Optional: Sort
    result.sort_by_score(false); // ascending for p-value

    result
}

/// Calculates enrichment statistics based on raw indices.
///
/// This function is designed for testing and scenarios where the ontology
/// structure is abstracted away into vectors of indices.
///
/// # Arguments
/// * `population_size` (N) - Total number of items in the population.
/// * `study_total` (n) - Total size of the study.
/// * `study_indices` - Indices of items in the study.
/// * `population_terms` - A slice of terms, where each term is a vector of item indices.
#[allow(non_snake_case)]
pub fn calculate_enrichment_from_indices(
    population_size: usize,
    study_total: usize,
    study_indices: &[usize],
    population_terms: &[Vec<usize>],
) -> Vec<PValue> {
    let cache = LogFactorialCache::new(population_size);
    let test = OneSidedEnrichmentTest;

    // Convert slice to HashSet for O(1) lookups during intersection
    let study_set: HashSet<usize> = study_indices.iter().cloned().collect();

    let mut measures: Vec<PValue> = Vec::with_capacity(population_terms.len());

    for term_indices in population_terms {
        let K = term_indices.len();

        // Calculate intersection (k)
        let k = term_indices
            .iter()
            .filter(|idx| study_set.contains(idx))
            .count();

        let hypergeo = Hypergeometric::new(study_total, K, population_size, &cache);
        let raw_p_value = test.calculate(&hypergeo, k);

        // Counts: (k=Annotated_Study, n=Total_Study, K=Annotated_Pop, N=Total_Pop)
        let counts = (
            k as u32,
            study_total as u32,
            K as u32,
            population_size as u32,
        );
        let measure = PValue::new(raw_p_value, counts);

        measures.push(measure);
    }

    // Apply Correction
    let correction = correction::None;
    correction.adjust(&mut measures, |m| &mut m.pvalue);

    measures
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::GeneSet;
    use crate::core::result::Measure;
    use csv::Writer;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use std::process;

    #[test]
    fn test_vector_based_logic() {
        let N = 10; // population gene count
        let n = 2; // gene gene count
        let study_indices = vec![1, 2];
        let terms = vec![
            vec![1, 2, 3], // Term A: K=3, k=2 (1,2 overlap)
            vec![1, 2],    // Term B: K=2, k=2 (1,2 overlap)
            vec![4, 5],    // Term C: K=2, k=0 (no overlap)
        ];

        let results = calculate_enrichment_from_indices(N, n, &study_indices, &terms);

        // Term A: k=2, n=2, K=3, N=10.
        // Hypergeometric(x=2, n=2, K=3, N=10)
        // Probability of drawing 2 specific items from 10 when target group is size 3.
        for result in &results {
            println!("P-value: {}", result.pvalue);
        }
        // Check counts structure
        // (k, n, K, N)
        assert_eq!(
            results[0].diagnostics().unwrap(),
            format!("{:?}", (2, 2, 3, 10))
        ); // Term A

        assert_eq!(
            results[1].diagnostics().unwrap(),
            format!("{:?}", (2, 2, 2, 10))
        ); // Term B

        assert_eq!(
            results[2].diagnostics().unwrap(),
            format!("{:?}", (0, 2, 2, 10))
        ); // Term C
    }

    #[test]
    fn test_go_analysis() {
        // Paths to the necessary files
        // These paths should be adjusted based on your local setup
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let population_genes_path = "tests/data/GOnone/population.txt";
        let study_genes_path = "tests/data/GOnone/study.txt";

        // ------ Load Gene Sets ------

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

        let study_genes =
            GeneSet::from_file(study_genes_path, &annotations).unwrap_or_else(|err| {
                eprintln!("Error: {}", err);
                process::exit(1);
            });

        let population_genes = GeneSet::from_file(population_genes_path, &annotations)
            .unwrap_or_else(|err| {
                eprintln!("Error: {}", err);
                process::exit(1);
            });

        let annotation_index = AnnotationIndex::new(
            annotations,
            &ontology,
            Some(&population_genes.recognized_genes()),
        );

        let result = analysis(&ontology, annotation_index, study_genes.recognized_genes());

        // Serialize to CSV
        let mut wtr = Writer::from_path("tests/data/GOnone/results_f.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }
}
