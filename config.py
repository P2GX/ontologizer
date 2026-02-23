import json

# Configuration matching the Problem struct in main.rs
config = {
    "method": "bayesian",
    "study_genes_path": "tests/data/GO0090717/study.txt",
    "population_genes_path": "tests/data/GO0090717/population.txt",
    "ontology_path": "tests/data/GO/go-basic.json",
    "annotation_path": "tests/data/GO/goa_human.gaf"
}

output_file = "problem.json"

with open(output_file, "w") as f:
    json.dump(config, f, indent=2)

print(f"{output_file} created successfully.")
