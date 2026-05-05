//! End-to-end frequentist enrichment on the yeast benchmark dataset.
//!
//! Runs the full pipeline (load ontology + annotations, build the index,
//! run the hypergeometric test with Benjamini-Hochberg correction) on the
//! yeast study/population sets shipped under `data/yeast/` and checks that
//! a meaningful fraction of the ground-truth GO terms in
//! `data/yeast/solution_yeast.tsv` lands in the top of the result list.
//!
//! Marked `#[ignore]` because it touches ~25 MB of GO/GAF files and takes a
//! few seconds. Run it explicitly with:
//!
//!     cargo test --release --test end_to_end_yeast -- --ignored --nocapture

use std::collections::HashSet;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use oboannotation::go::{GoAnnotations, GoGafAnnotationLoader};
use oboannotation::io::AnnotationLoader;
use ontolius::io::OntologyLoaderBuilder;
use ontolius::ontology::csr::FullCsrOntology;

use ontologizer::{AnnotationIndex, Background, Correction, GeneSet, frequentist_analysis};

const GO_URL: &str = "https://purl.obolibrary.org/obo/go/go-basic.json";
const GAF_URL: &str = "https://current.geneontology.org/annotations/sgd.gaf.gz";

/// Top-K result slice we examine for ground-truth recovery.
const TOP_K: usize = 50;

/// Minimum number of solution terms we expect in the top-K. Calibrated against
/// the BH-corrected frequentist run on the yeast dataset; leave a comfortable
/// margin so ontology updates don't cause spurious failures.
const MIN_RECALL: usize = 5;

fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn ensure_cached(url: &str, target: &Path, gzipped: bool) -> io::Result<()> {
    if target.exists() {
        return Ok(());
    }
    if let Some(parent) = target.parent() {
        fs::create_dir_all(parent)?;
    }
    let response = ureq::get(url)
        .call()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    let body = response.into_body().into_reader();
    let mut reader: Box<dyn io::Read> = if gzipped {
        Box::new(GzDecoder::new(body))
    } else {
        Box::new(body)
    };
    let mut writer = io::BufWriter::new(fs::File::create(target)?);
    io::copy(&mut reader, &mut writer)?;
    Ok(())
}

fn read_solution_terms(path: &Path) -> HashSet<String> {
    let raw = fs::read_to_string(path).expect("solution file missing");
    raw.lines()
        .filter(|l| !l.trim().is_empty())
        .filter_map(|l| l.split('\t').next().map(|s| s.trim().to_string()))
        .collect()
}

#[test]
#[ignore = "Downloads ~25 MB of ontology data and runs the full pipeline"]
fn yeast_frequentist_recovers_truth() {
    let root = project_root();
    let go_path = root.join("go-basic.json");
    let gaf_path = root.join("goa_yeast.gaf");
    let study_path = root.join("data/yeast/study_genes_yeast.txt");
    let pop_path = root.join("data/yeast/population_genes_yeast.txt");
    let solution_path = root.join("data/yeast/solution_yeast.tsv");

    ensure_cached(GO_URL, &go_path, false).expect("download go-basic.json");
    ensure_cached(GAF_URL, &gaf_path, true).expect("download yeast GAF");

    let ontology: FullCsrOntology = OntologyLoaderBuilder::new()
        .obographs_parser()
        .build()
        .load_from_path(&go_path)
        .expect("load ontology");
    let annotations: GoAnnotations = GoGafAnnotationLoader
        .load_from_path(&gaf_path)
        .expect("load GAF");

    let population = GeneSet::from_file(pop_path.to_str().unwrap(), None).expect("load population");
    let study =
        GeneSet::from_file(study_path.to_str().unwrap(), Some(&population)).expect("load study");

    let index = AnnotationIndex::new(annotations, &ontology, population.recognized_genes());

    let result = frequentist_analysis(
        &ontology,
        &index,
        study.recognized_genes(),
        &Background::Standard,
        &Correction::BenjaminiHochberg,
    );

    let truth = read_solution_terms(&solution_path);
    assert!(!truth.is_empty(), "solution file produced no terms");

    let top: HashSet<&str> = result
        .iter_items()
        .take(TOP_K)
        .map(|item| item.id.as_str())
        .collect();
    let hits: Vec<&str> = top
        .iter()
        .copied()
        .filter(|id| truth.contains(*id))
        .collect();

    println!(
        "yeast frequentist + BH: {} / {} truth terms recovered in top {}",
        hits.len(),
        truth.len(),
        TOP_K
    );
    println!("hits: {:?}", hits);

    assert!(
        hits.len() >= MIN_RECALL,
        "expected ≥{MIN_RECALL} of {} truth terms in top {TOP_K}, got {}",
        truth.len(),
        hits.len(),
    );
}
