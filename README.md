# Ontologizer

Fast and safe implementation of the Ontologizer — a tool for Gene Ontology (GO) enrichment analysis using Frequentist (
hypergeometric test) and Bayesian (inference) methods.

## Gene Symbols

Gene symbols in the study and population gene sets must match the gene symbols in the `.gaf` annotation file —
specifically column 2 (`DB_Object_Symbol`).

## Try it

The `examples/` folder ships ready-to-run study/population gene sets for **human** and **yeast**, with two configs per
organism — one Frequentist, one Bayesian. The only file you need to provide yourself is the organism-specific GO
annotation (`.gaf`); the GO ontology JSON is fetched automatically on first run.

### 1. Build

```bash
cargo build --release
```

### 2. Download the GAF for your organism

Place the unzipped file in `examples/` with the filename the configs expect:

**Human** (`examples/goa_human.gaf`):

```bash
wget https://current.geneontology.org/annotations/goa_human.gaf.gz
gunzip goa_human.gaf.gz
mv goa_human.gaf examples/goa_human.gaf
```

**Yeast** (`examples/goa_yeast.gaf` — the SGD-curated yeast GAF, renamed to match the config):

```bash
wget https://current.geneontology.org/annotations/sgd.gaf.gz
gunzip sgd.gaf.gz
mv sgd.gaf examples/goa_yeast.gaf
```

### 3. Run

Pick one (or several) of the four configs:

```bash
cargo run --release -- examples/human/config_frequentist.json
cargo run --release -- examples/human/config_bayesian.json
cargo run --release -- examples/yeast/config_frequentist.json
cargo run --release -- examples/yeast/config_bayesian.json
```

Each invocation writes its results next to the inputs, e.g.
`examples/human/results_human_frequentist.csv`,
`examples/yeast/results_yeast_bayesian.csv`.

### 4. Sanity-check (optional)

Each example folder includes a `solution_{org}.tsv` file mapping known-enriched GO terms to their gene members. It is a
ground-truth reference, not a results CSV — useful for spot-checking that top hits in your output overlap with these
terms.

## Configuration

Each example comes with a ready-to-use `config.json`. The schema is:

```json
{
  "study_file": "examples/ORGANISM/NAME_study_genes.txt",
  "pop_file": "examples/ORGANISM/NAME_pop_genes.txt",
  "go_file": "examples/go-basic.json",
  "goa_file": "examples/goa_yeast.gaf",
  "out_file": "examples/ORGANISM/results_ORGANISM_bayesian.csv",
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