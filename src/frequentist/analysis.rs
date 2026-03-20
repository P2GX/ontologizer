use indexmap::IndexSet;
use std::collections::HashSet;

use crate::core::AnnotationIndex;
use crate::core::result::AnalysisResult;
use crate::frequentist;
use crate::frequentist::algorithm::{OneSidedEnrichmentTest, StatisticalTest};
use crate::frequentist::background::Restrict;
use crate::frequentist::correction::Adjustment;
use crate::frequentist::distribution::{Hypergeometric, LogFactorialCache};
use crate::frequentist::measure::PValue;
use ontolius::ontology::csr::FullCsrOntology;

#[allow(non_snake_case)]
pub fn analysis(
    ontology: &FullCsrOntology,
    annotations: &AnnotationIndex,
    study_genes: &HashSet<String>,
    topology: &frequentist::Background,
    correction: &frequentist::Correction,
) -> AnalysisResult {
    let N_pop = annotations.get_genes().len();
    let cache = LogFactorialCache::new(N_pop);
    let test = OneSidedEnrichmentTest;

    let study_indices: IndexSet<usize> = study_genes
        .iter()
        .filter_map(|gene| annotations.get_idx_by_gene(gene))
        .collect();

    let mut measures: Vec<PValue> = Vec::new();

    for (term_idx, _) in annotations.get_terms().iter().enumerate() {
        let parent_indices = topology.restrict(term_idx, annotations, ontology);
        let n = study_indices.intersection(&parent_indices).count();
        let N = parent_indices.len();

        let genes_for_term = annotations.get_gene_idxs_for_term_idx(term_idx);
        let k = study_indices.intersection(genes_for_term).count();
        let K = genes_for_term.intersection(&parent_indices).count();

        let hypergeo = Hypergeometric::new(n, K, N, &cache);
        let raw_p_value = test.calculate(&hypergeo, k);

        let counts = (k as u32, n as u32, K as u32, N as u32);
        measures.push(PValue::new(raw_p_value, counts));
    }

    correction.adjust(&mut measures, |m| &mut m.pvalue);

    let observed_genes: Vec<bool> = (0..N_pop).map(|i| study_indices.contains(&i)).collect();

    let mut result =
        AnalysisResult::from_measures(&measures, ontology, annotations, &observed_genes);
    result.sort_by_score(false);

    let topology = format!("{:?}", topology);
    let correction = format!("{:?}", correction);

    result.with_meta(&[
        ("Method", "Frequentist"),
        ("Topology", &topology),
        ("Correction", &correction),
    ])
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Background;
    use crate::Correction;
    use crate::GeneSet;
    use crate::core::result::Measure;
    use csv::Writer;
    use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
    use oboannotation::io::AnnotationLoader;
    use ontolius::io::OntologyLoaderBuilder;
    use std::process;

    #[allow(non_snake_case)]
    fn calculate_enrichment_from_indices(
        population_size: usize,
        study_total: usize,
        study_indices: &[usize],
        population_terms: &[Vec<usize>],
        correction: &Correction,
    ) -> Vec<PValue> {
        let cache = LogFactorialCache::new(population_size);
        let test = OneSidedEnrichmentTest;

        let study_set: HashSet<usize> = study_indices.iter().cloned().collect();

        let mut measures: Vec<PValue> = Vec::with_capacity(population_terms.len());

        for term_indices in population_terms {
            let K = term_indices.len();
            let k = term_indices
                .iter()
                .filter(|idx| study_set.contains(idx))
                .count();

            let hypergeo = Hypergeometric::new(study_total, K, population_size, &cache);
            let raw_p_value = test.calculate(&hypergeo, k);

            let counts = (
                k as u32,
                study_total as u32,
                K as u32,
                population_size as u32,
            );
            measures.push(PValue::new(raw_p_value, counts));
        }

        correction.adjust(&mut measures, |m| &mut m.pvalue);

        measures
    }

    #[test]
    fn test_vector_based_logic() {
        #[allow(non_snake_case)]
        let N = 10;
        let n = 2;
        let study_indices = vec![1, 2];
        let terms = vec![
            vec![1, 2, 3], // Term A: K=3, k=2 (1,2 overlap)
            vec![1, 2],    // Term B: K=2, k=2 (1,2 overlap)
            vec![4, 5],    // Term C: K=2, k=0 (no overlap)
        ];

        let results =
            calculate_enrichment_from_indices(N, n, &study_indices, &terms, &Correction::None);

        for result in &results {
            println!("P-value: {}", result.pvalue);
        }

        assert_eq!(
            results[0].diagnostics().unwrap(),
            format!("{:?}", (2, 2, 3, 10))
        );
        assert_eq!(
            results[1].diagnostics().unwrap(),
            format!("{:?}", (2, 2, 2, 10))
        );
        assert_eq!(
            results[2].diagnostics().unwrap(),
            format!("{:?}", (0, 2, 2, 10))
        );
    }

    #[test]
    fn test_go_analysis() {
        let go_path = "tests/data/GO/go-basic.json";
        let gaf_path = "tests/data/GO/goa_human.gaf";
        let population_genes_path = "tests/data/GOnone/population.txt";
        let study_genes_path = "tests/data/GOnone/study.txt";

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

        let result = analysis(
            &ontology,
            &annotation_index,
            &study_genes.recognized_genes(),
            &Background::Standard,
            &Correction::None,
        );

        let mut wtr = Writer::from_path("tests/data/GOnone/results_f.csv").unwrap();
        for item in result.items {
            wtr.serialize(item).unwrap();
        }
    }
}
