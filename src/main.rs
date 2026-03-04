use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;
use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;
use ontologizer::{AnnotationIndex, EnrichmentResult, GeneSet};

use csv::Writer;
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process;

// Domain-specific URLs for the Gene Ontology (GO) and Annotations.
const GO_URL: &str = "https://purl.obolibrary.org/obo/go/go-basic.json";
const GOA_BASE_URL: &str = "http://current.geneontology.org/annotations/";

fn map_organism_to_gaf(organism: &str) -> String {
    format!("{}goa_{}.gaf.gz", GOA_BASE_URL, organism.to_lowercase())
}

fn ensure_cached(
    url: &str,
    target_path: &Path,
    is_gzipped: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if target_path.exists() {
        println!("Using existing file: {:?}", target_path);
        return Ok(());
    }

    println!("Downloading {} to {:?}", url, target_path);

    if let Some(parent_dir) = target_path.parent() {
        fs::create_dir_all(parent_dir)?;
    }

    let response: ureq::http::Response<ureq::Body> = ureq::get(url).call()?;
    let mut reader: Box<dyn io::Read> = if is_gzipped {
        Box::new(GzDecoder::new(response.into_body().into_reader()))
    } else {
        Box::new(response.into_body().into_reader())
    };

    let file: File = File::create(target_path)?;
    let mut writer: BufWriter<File> = BufWriter::new(file);

    io::copy(&mut reader, &mut writer)?;

    Ok(())
}

#[derive(Deserialize, Debug)]
pub struct Config {
    pub organism: String,
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
        println!("No config file provided. Loading from 'problem.json'");
        "problem.json"
    };

    // Open the configuration file
    let file = File::open(config_path).unwrap_or_else(|err| {
        eprintln!("Error opening config file '{}': {}", config_path, err);
        eprintln!("Usage: cargo run -- [path/to/problem.json]");
        process::exit(1);
    });
    let reader = BufReader::new(file);

    let config: Config = serde_json::from_reader(reader).unwrap_or_else(|err| {
        eprintln!("Error parsing JSON config: {}", err);
        process::exit(1);
    });

    // ------ Local Data Directory Resolution ------
    let current_dir: PathBuf = env::current_dir().unwrap_or_else(|err| {
        eprintln!("Failed to determine current directory: {}", err);
        process::exit(1);
    });
    let data_dir: PathBuf = current_dir.join("data");
    fs::create_dir_all(&data_dir).unwrap_or_else(|err| {
        eprintln!("Failed to create data directory at {:?}: {}", data_dir, err);
        process::exit(1);
    });

    // ------ Default Data Locations ------
    let ontology_path: PathBuf = data_dir.join("go-basic.json");
    let annotation_path: PathBuf = data_dir.join(format!("goa_{}.gaf", config.organism));

    // ------ Ensure Gene Ontology and Annotations exist ------
    ensure_cached(GO_URL, &ontology_path, false).unwrap_or_else(|err| {
        eprintln!("Failed to fetch ontology: {}", err);
        process::exit(1);
    });
    let goa_url: String = map_organism_to_gaf(&config.organism);
    ensure_cached(&goa_url, &annotation_path, true).unwrap_or_else(|err| {
        eprintln!("Failed to fetch annotation: {}", err);
        process::exit(1);
    });

    // ------ Load Gene Ontology and Annotations ------
    println!("Processing Ontology and Annotation files...");
    let ontology_loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let ontology: FullCsrOntology = ontology_loader
        .load_from_path(&ontology_path)
        .unwrap_or_else(|err| {
            eprintln!("Error loading ontology: {}", err);
            process::exit(1);
        });

    let annotations_loader = GoGafAnnotationLoader;
    let annotations: GoAnnotations = annotations_loader
        .load_from_path(&annotation_path)
        .unwrap_or_else(|err| {
            eprintln!(
                "Could not load GAF file from {:?}: {}",
                annotation_path, err
            );
            process::exit(1);
        });

    // ------ Load Gene Sets ------
    println!("Loading study gene set from {:?}", config.study_genes_path);
    let study_genes =
        GeneSet::from_file(&config.study_genes_path, &annotations).unwrap_or_else(|err| {
            eprintln!("Error loading study genes: {}", err);
            process::exit(1);
        });

    println!(
        "Loading population gene set from {:?}",
        config.population_genes_path
    );
    let population_genes = GeneSet::from_file(&config.population_genes_path, &annotations)
        .unwrap_or_else(|err| {
            eprintln!("Error loading population genes: {}", err);
            process::exit(1);
        });

    let annotation_index = AnnotationIndex::new(
        annotations,
        &ontology,
        Some(&population_genes.recognized_genes()),
    );

    println!("Starting enrichment analysis...");
    let result: EnrichmentResult = match config.method {
        MethodConfig::Frequentist => ontologizer::frequentist_analysis(
            &ontology,
            &annotation_index,
            &study_genes.recognized_genes(),
        ),
        MethodConfig::Bayesian => ontologizer::bayesian_analysis(
            &ontology,
            &annotation_index,
            &study_genes.recognized_genes(),
        ),
    };

    // Serialize to CSV
    let result_dir: PathBuf = current_dir.join("result");
    fs::create_dir_all(&result_dir).unwrap_or_else(|err| {
        eprintln!(
            "Failed to create output directory at {:?}: {}",
            result_dir, err
        );
        process::exit(1);
    });

    let output_filename: PathBuf = result_dir.join(format!("{:?}_results.csv", config.method));
    println!("Writing results to: {:?}", output_filename);

    let mut wtr = Writer::from_path(&output_filename).unwrap_or_else(|err| {
        eprintln!("Error creating output file {:?}: {}", output_filename, err);
        process::exit(1);
    });

    for item in result.items {
        wtr.serialize(item).unwrap_or_else(|err| {
            eprintln!("Error serializing record: {}", err);
            process::exit(1);
        });
    }

    println!("Done!");
}
