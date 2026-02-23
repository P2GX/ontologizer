use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;

use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;

use csv::Writer;
use ontologizer;
use ontologizer::{AnnotationIndex, EnrichmentResult, GeneSet};
use serde::Deserialize;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::process;

#[derive(Deserialize, Debug)]
pub struct Config {
    pub ontology_path: String,
    pub annotation_path: String,
    pub study_genes_path: String,
    pub population_genes_path: String,
    #[serde(flatten)]
    pub method: MethodConfig,
}

#[derive(Deserialize, Debug)]
#[serde(tag = "method", rename_all = "lowercase")]
pub enum MethodConfig {
    Frequentist,
    Bayesian,
}
pub fn main() {
    let args: Vec<String> = env::args().collect();

    let config_path = if args.len() > 1 {
        &args[1]
    } else {
        println!("No config file provided. Loading from '/problem.json'");
        "problem.json"
    };

    // Open the configuration file
    let file = File::open(config_path).unwrap_or_else(|err| {
        eprintln!("Error opening config file '{}': {}", config_path, err);
        eprintln!("Usage: cargo run -- [path/to/problem.json]");
        process::exit(1);
    });
    let reader = BufReader::new(file);

    let problem: Config = serde_json::from_reader(reader).unwrap_or_else(|err| {
        eprintln!("Error parsing JSON config: {}", err);
        process::exit(1);
    });

    // ------ Load Gene Ontology and Annotations ------
    let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let ontology: FullCsrOntology = ontology_loader
        .load_from_path(&problem.ontology_path)
        .unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

    let annotations_loader = GoGafAnnotationLoader;
    let annotations: GoAnnotations = annotations_loader
        .load_from_path(&problem.annotation_path)
        .expect("Could not load GAF file");

    // ------ Load Gene Sets ------
    let study_genes =
        GeneSet::from_file(&problem.study_genes_path, &annotations).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });
    let population_genes = GeneSet::from_file(&problem.population_genes_path, &annotations)
        .unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

    let annotation_index = AnnotationIndex::new(
        annotations,
        &ontology,
        Some(&population_genes.recognized_genes()),
    );

    let result: EnrichmentResult = match problem.method {
        MethodConfig::Frequentist => ontologizer::frequentist_analysis(
            &ontology,
            annotation_index,
            study_genes.recognized_genes(),
        ),
        MethodConfig::Bayesian => ontologizer::bayesian_analysis(
            &ontology,
            annotation_index,
            study_genes.recognized_genes(),
        ),
    };

    // Serialize to CSV
    let output_filename = format!("result/{:?}_results.csv", problem.method);
    let mut wtr = Writer::from_path(&output_filename).unwrap();
    for item in result.items {
        wtr.serialize(item).unwrap();
    }

    println!("Done! Results written to {}", output_filename);
}
