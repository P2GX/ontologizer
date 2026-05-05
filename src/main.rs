use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;
use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;
use ontologizer::{AnalysisResult, AnnotationIndex, GeneSet, Method};

use flate2::read::GzDecoder;
use serde::Deserialize;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process;

// Domain-specific URLs for the Gene Ontology (GO) and Annotations.
const GO_URL: &str = "https://purl.obolibrary.org/obo/go/go-basic.json";

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
    pub study_file: String,
    pub pop_file: String,
    pub go_file: String,
    pub goa_file: String,
    pub out_file: String,
    pub method: Method,
}

pub fn main() {
    let args: Vec<String> = env::args().collect();

    let config_file = if args.len() > 1 {
        &args[1]
    } else {
        eprintln!("Error: no config file provided.");
        eprintln!("Usage: cargo run -- path/to/config.json");
        process::exit(1);
    };

    // Open the configuration file
    let file = File::open(config_file).unwrap_or_else(|err| {
        eprintln!("Error opening config file '{}': {}", config_file, err);
        eprintln!("Usage: cargo run -- [path/to/config.json]");
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

    let out_file: PathBuf = PathBuf::from(&config.out_file);
    let out_dir = out_file.parent().unwrap_or_else(|| {
        eprintln!(
            "No parent directoy for Output file at: {:?}.",
            &config.out_file
        );
        process::exit(1);
    });

    fs::create_dir_all(out_dir).unwrap_or_else(|err| {
        eprintln!(
            "Failed to create output directory at {:?}: {}",
            &config.out_file, err
        );
        process::exit(1);
    });

    // ------ Default Data Locations ------
    let ontology_path: PathBuf = current_dir.join(&config.go_file);
    let annotation_path: PathBuf = current_dir.join(&config.goa_file);

    // ------ Ensure Gene Ontology exist ------
    ensure_cached(GO_URL, &ontology_path, false).unwrap_or_else(|err| {
        eprintln!("Failed to fetch ontology: {}", err);
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
    println!("Loading population gene set from {:?}", config.pop_file);
    let population_genes = GeneSet::from_file(&config.pop_file, None).unwrap_or_else(|err| {
        eprintln!("Error loading population genes: {}", err);
        process::exit(1);
    });

    println!("Loading study gene set from {:?}", config.study_file);
    let study_genes = GeneSet::from_file(&config.study_file, Some(&population_genes))
        .unwrap_or_else(|err| {
            eprintln!("Error loading study genes: {}", err);
            process::exit(1);
        });

    let annotation_index =
        AnnotationIndex::new(annotations, &ontology, population_genes.recognized_genes());

    println!("Starting enrichment analysis...");
    let results: AnalysisResult = match config.method {
        Method::Frequentist {
            background,
            correction,
        } => ontologizer::frequentist_analysis(
            &ontology,
            &annotation_index,
            &study_genes.recognized_genes(),
            &background,
            &correction,
        ),
        Method::Bayesian => ontologizer::bayesian_analysis(
            &ontology,
            &annotation_index,
            &study_genes.recognized_genes(),
        ),
    };

    // Serialize to CSV
    println!("Writing results to: {:?}", out_file);

    results.save_to_csv(&out_file, true).unwrap_or_else(|err| {
        eprintln!("Error writing output file {:?}: {}", out_file, err);
        process::exit(1);
    });

    println!("Done!");
}
