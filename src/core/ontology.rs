use ontolius::{
    io::OntologyLoaderBuilder,
    ontology::{MetadataAware, OntologyTerms, csr::FullCsrOntology},
};

pub struct Ontology {
    ontology: FullCsrOntology,
}

impl Ontology {
    pub fn new(go_json: &str) -> Self {
        let loader = OntologyLoaderBuilder::new().obographs_parser().build();
        let go: FullCsrOntology = loader
            .load_from_path(go_json)
            .expect("Could not load ontology");
        Self { ontology: go }
    }

    pub fn term_count(&self) -> usize {
        self.ontology.len()
    }

    pub fn ontology(&self) -> &FullCsrOntology {
        &self.ontology
    }

    pub fn version(&self) -> String {
        self.ontology.version().to_string()
    }
}
