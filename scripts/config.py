import json
import pathlib
from typing import Any, Dict, List, Optional

def generate_config(
    study_genes_path: str,
    population_genes_path: str,
    go_path: str,
    goa_path: str,
    method: str,
    topology: Optional[str] = None,
    correction: Optional[str] = None,
) -> None:
    """
    Generates a nested JSON configuration for the Rust Config struct.
    Target output is ../config.json relative to this script's location.
    """

    # Resolve output path: script in /project/scripts/, target is /project/config.json
    script_dir: pathlib.Path = pathlib.Path(__file__).parent.resolve()
    output_path: pathlib.Path = script_dir.parent / "config.json"

    # Base configuration mapping
    config: Dict[str, Any] = {
        "study_genes_path": study_genes_path,
        "population_genes_path": population_genes_path,
        "go_path": go_path,
        "goa_path": goa_path,
    }

    # Method logic: Nested construction to match Rust Enum (non-flattened)
    method_lower: str = method.lower()

    if method_lower == "frequentist":
        valid_topologies: List[str] = ["Standard", "ParentUnion", "ParentIntersection"]
        valid_corrections: List[str] = ["Bonferroni", "BonferroniHolm", "BenjaminHochberg", "None"]

        if topology not in valid_topologies:
            raise ValueError(f"Invalid topology '{topology}'. Must be one of: {valid_topologies}")

        if correction not in valid_corrections:
            raise ValueError(f"Invalid correction '{correction}'. Must be one of: {valid_corrections}")

        config["method"] = {
            "method": "frequentist",
            "topology": topology,
            "correction": correction
        }
    elif method_lower == "bayesian":
        config["method"] = {
            "method": "bayesian"
        }
    else:
        raise ValueError("Method must be 'frequentist' or 'bayesian'")

    # Atomic write to file
    with open(output_path, "w") as f:
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    # Example execution within the script
    generate_config(
        study_genes_path="data/Genes_Yeast/study_genes.txt",
        population_genes_path="data/Genes_Yeast/population_genes.txt",
        go_path="data/GO/go-basic.json",
        goa_path="data/GOA/gene_association.sgd.20260120.gaf",
        method="bayesian",
        topology="Standard",
        correction="None"
    )