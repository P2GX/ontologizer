use ontolius::{
    io::OntologyLoaderBuilder,
    ontology::{OntologyTerms, csr::FullCsrOntology},
};

pub struct Ontologizer {
    ontology: FullCsrOntology,
}

impl Ontologizer {
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
}
// region:    --- Tests

#[cfg(test)]
mod tests {
    type Result<T> = core::result::Result<T, Box<dyn std::error::Error>>; // For tests.

    use super::*;

    #[test]
    #[ignore]
    fn test_name() -> Result<()> {
        let go_path = "/Users/robin/data/go.json";
        let ontologizer = Ontologizer::new(go_path);

        assert!(ontologizer.term_count() == 7);

        Ok(())
    }

    #[test]
    #[ignore]
    fn test2() -> Result<()> {
        let x = 2;
        let y = 3;

        assert!(x + y == 5);

        Ok(())
    }
}

// endregion: --- Tests
