use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;

use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;

use csv::Writer;
use serde::Deserialize;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::process;

mod bayesian;
mod core;
mod frequentist;

#[derive(Deserialize, Debug)]
struct Problem {
    method: String,
    study_genes_path: String,
    population_genes_path: String,
    ontology_path: String,
    annotation_path: String,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let config_path = if args.len() > 1 {
        &args[1]
    } else {
        "problem.json"
    };

    // Open the configuration file
    let file = File::open(config_path).unwrap_or_else(|err| {
        eprintln!("Error opening config file '{}': {}", config_path, err);
        eprintln!("Usage: cargo run -- [path/to/problem.json]");
        process::exit(1);
    });
    let reader = BufReader::new(file);

    let problem: Problem = serde_json::from_reader(reader).unwrap_or_else(|err| {
        eprintln!("Error parsing JSON config: {}", err);
        process::exit(1);
    });

    if problem.method != "frequentist" && problem.method != "bayesian" {
        eprintln!(
            "Error: Method argument must be `frequentist` or `bayesian`, got `{}`.",
            problem.method
        );
        process::exit(1);
    }

    // ------ Load Gene Sets ------
    let study_genes = core::load_gene_set(&problem.study_genes_path).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    let population_genes =
        core::load_gene_set(&problem.population_genes_path).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
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

    let annotation_index =
        core::AnnotationIndex::new(annotations, &ontology, Some(&population_genes));

    let result;
    if problem.method == "frequentist" {
        result = frequentist::run(&ontology, annotation_index, study_genes);
    } else {
        result = bayesian::run(&ontology, annotation_index, study_genes)
    }

    // Serialize to CSV
    let output_filename = format!("result/{}_results.csv", problem.method);
    let mut wtr = Writer::from_path(&output_filename).unwrap();
    for item in result.items {
        wtr.serialize(item).unwrap();
    }

    println!("Done! Results written to {}", output_filename);
}
