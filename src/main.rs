use ontolius::ontology::csr::FullCsrOntology;
use ontolius::io::OntologyLoaderBuilder;

use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;

use std::env;
use std::process;
use csv::Writer;

mod core;
mod bayesian;
mod frequentist;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        eprintln!("Error: Expected 5 arguments, got {}.", args.len() - 1);
        eprintln!("Usage: cargo run -- <method> <path/to/study-gene-set> <path/to/population-gene-set> <path/to/gene-ontology> <path/to/gene-ontology-annotations");
        process::exit(1);
    }

    let method = &args[1];
    if method != "frequentist" && method != "bayesian" {
        eprintln!(
            "Error: Method argument must be `frequentist` or `bayesian`, got `{}`.",
            method
        );
        process::exit(1);
    }

    // ------ Load Gene Sets ------
    let study_genes_path = &args[2];
    let study_genes = core::load_gene_set(study_genes_path).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    let population_genes_path = &args[3];
    let population_genes = core::load_gene_set(population_genes_path).unwrap_or_else(|err| {
        eprintln!("Error: {}", err);
        process::exit(1);
    });

    // ------ Load Gene Ontology and Annotations ------
    let go_path = &args[4];
    let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let ontology: FullCsrOntology = ontology_loader
        .load_from_path(go_path)
        .unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            process::exit(1);
        });

    let gaf_path = &args[5];
    let annotations_loader = GoGafAnnotationLoader;
    let annotations: GoAnnotations = annotations_loader
        .load_from_path(gaf_path)
        .expect("Could not load GAF file");
    
    
    let annotation_index = core::AnnotationIndex::new(annotations, &ontology, Some(&population_genes));

    let result;
    if method == "frequentist"{
        result = bayesian::run(&ontology, annotation_index, study_genes);
    }
    else {
        result = bayesian::run(&ontology, annotation_index, study_genes)
    }


    // Serialize to CSV
    let mut wtr = Writer::from_path("results.csv").unwrap();
    for item in result.items {
        wtr.serialize(item).unwrap();
    }

    dbg!(args);
}
