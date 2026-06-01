# Ontologizer

Fast and safe implementation of the Ontologizer — a tool for Gene Ontology (GO) enrichment analysis using Frequentist
(hypergeometric test) or Bayesian (inference) methods.

## Gene Symbols

Gene symbols in the study and population gene sets must match the gene symbols in the `.gaf` annotation file —
specifically column 2 (`DB_Object_Symbol`).

## Try it

The `examples/` folder ships ready-to-run study/population gene sets for **yeast**, **human**, **fly**, **mouse**, and
**rat**. The GO ontology and the per-organism GAF must be downloaded manually (see step 2).

### 1. Build

```bash
cargo build --release
```

### 2. Download the GO ontology and the GAF for your organism

Drop the GO ontology and the gzipped GAF into `examples/` under the filenames the config expects — the binary reads
`.gaf.gz` directly, no unzip step needed.

```bash
wget -P examples https://purl.obolibrary.org/obo/go/go-basic.json
```

| Organism | Target path                 | Source                                                          |
|----------|-----------------------------|-----------------------------------------------------------------|
| yeast    | `examples/sgd.gaf.gz`       | `https://current.geneontology.org/annotations/sgd.gaf.gz`       |
| human    | `examples/goa_human.gaf.gz` | `https://current.geneontology.org/annotations/goa_human.gaf.gz` |
| fly      | `examples/fb.gaf.gz`        | `https://current.geneontology.org/annotations/fb.gaf.gz`        |
| mouse    | `examples/mgi.gaf.gz`       | `https://current.geneontology.org/annotations/mgi.gaf.gz`       |
| rat      | `examples/rgd.gaf.gz`       | `https://current.geneontology.org/annotations/rgd.gaf.gz`       |

For example:

```bash
wget -P examples https://current.geneontology.org/annotations/goa_human.gaf.gz
```

### 3. Generate a config

```bash
python3 examples/generate_config.py <organism> <method> [--correction <correction>]
```

For example:

```bash
python3 examples/generate_config.py human bayesian
python3 examples/generate_config.py yeast frequentist --correction Bonferroni
```

This writes `examples/<organism>/config.json`.

### 4. Run

```bash
cargo run --release -- examples/<organism>/config.json
```

Results are written next to the inputs as `examples/<organism>/results_<organism>.tsv`.

### 5. Sanity-check (optional)

Each example folder includes a `solution_{org}.tsv` file mapping known-enriched GO terms to their gene members. It is a
ground-truth reference, not a results CSV — useful for spot-checking that top hits in your output overlap with these
terms.

## Configuration

`examples/generate_config.py` produces a `config.json` per organism. The schema is:

```json
{
  "study_file": "examples/ORGANISM/study_genes_ORGANISM.txt",
  "pop_file": "examples/ORGANISM/population_genes_ORGANISM.txt",
  "go_file": "examples/go-basic.json",
  "goa_file": "examples/GAF_FILENAME.gaf.gz",
  "out_file": "examples/ORGANISM/results_ORGANISM.tsv",
  "method": {
    "method": "Bayesian"
  }
}
```

**Method otions:** `Bayesian`, `Frequentist`.

For frequentist analysis a further parameter must be passed

```json
{
  "method": {
    "method": "Frequentist",
    "correction": "Bonferroni"
  }
}
```

**Correction options:** `Bonferroni`, `BonferroniHolm`, `BenjaminiHochberg`, `None`