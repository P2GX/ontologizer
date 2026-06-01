import os
import json
import argparse
from typing import Any, Dict, List, Optional


# Filename each organism's GAF is expected to live under, in examples/.
# Download instructions are in README.md.
GAF_FILE: Dict[str, str] = {
    "yeast": "sgd.gaf.gz",
    "human": "goa_human.gaf.gz",
    "fly":   "fb.gaf.gz",
    "mouse": "mgi.gaf.gz",
    "rat":   "rgd.gaf.gz",
}


def generate_config(
    study_file: str,
    pop_file: str,
    go_file: str,
    goa_file: str,
    out_file: str,
    method: str,
    output_path: str,
    correction: Optional[str] = None,
) -> None:
    """
    Generates a JSON configuration for the Rust `Config` struct in src/main.rs
    and writes it to `output_path`.
    """
    config: Dict[str, Any] = {
        "study_file": study_file,
        "pop_file": pop_file,
        "go_file": go_file,
        "goa_file": goa_file,
        "out_file": out_file,
    }

    # `Method` is an internally-tagged enum (`#[serde(tag = "method")]`),
    # so the variant name appears verbatim under the "method" key.
    method_lower: str = method.lower()

    if method_lower == "frequentist":
        valid_corrections: List[str] = [
            "Bonferroni",
            "BonferroniHolm",
            "BenjaminiHochberg",
            "None",
        ]
        if correction not in valid_corrections:
            raise ValueError(
                f"Invalid correction '{correction}'. Must be one of: {valid_corrections}"
            )
        config["method"] = {
            "method": "Frequentist",
            "correction": correction,
        }
    elif method_lower == "bayesian":
        config["method"] = {
            "method": "Bayesian",
        }
    else:
        raise ValueError("Method must be 'frequentist' or 'bayesian'")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(config, f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("organism", choices=list(GAF_FILE))
    parser.add_argument("method", choices=["bayesian", "frequentist"])
    parser.add_argument(
        "--mtc", "--correction",
        choices=["Bonferroni", "BonferroniHolm", "BenjaminiHochberg", "None"],
        default="Bonferroni",
    )
    args = parser.parse_args()

    org: str = args.organism
    gaf_name: str = GAF_FILE[org]

    repo_root: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_path: str = os.path.join(repo_root, "examples", org, "config.json")

    generate_config(
        study_file=f"examples/{org}/study_genes_{org}.txt",
        pop_file=f"examples/{org}/population_genes_{org}.txt",
        go_file="examples/go-basic.json",
        goa_file=f"examples/{gaf_name}",
        out_file=f"examples/{org}/results_{org}.tsv",
        method=args.method,
        output_path=output_path,
        correction=args.mtc if args.method == "frequentist" else None,
    )