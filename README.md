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

> **Heads up:** gene symbols in study/pop files must match column 2 (`DB_Object_Symbol`) of the GAF. Organism mismatch
> is the most common first-run failure.

### 1. Build

```bash
cargo build --release
```

### 2. Download the GAF for your organism

Place the unzipped file in `data/` with the filename the configs expect:

**Human** (`data/goa_human.gaf`):

```bash
mkdir -p data
wget https://current.geneontology.org/annotations/goa_human.gaf.gz
gunzip goa_human.gaf.gz
mv goa_human.gaf data/
```

**Yeast** (`data/goa_yeast.gaf` — the SGD-curated yeast GAF, renamed to match the config):

```bash
mkdir -p data
wget https://current.geneontology.org/annotations/sgd.gaf.gz
gunzip sgd.gaf.gz
mv sgd.gaf data/goa_yeast.gaf
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

Running with no argument falls back to the project-root `config.json` (yeast Frequentist):

```bash
cargo run --release
```

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
  "go_file": "data/go-basic.json",
  "goa_file": "data/goa_yeast.gaf",
  "out_file": "path/to/output.csv",
  "method": {
    "method": "Bayesian"
  }
}
```

**Method otions:** `Bayesian`, `Frequentist`.

For frequentist analysis two further parameters must be passed

```json
{
  "method": {
    "method": "Frequentist",
    "background": "Standard",
    "correction": "Bonferroni"
  }
}
```

**Background options:** `Standard`, `ParentUnion`, `ParentIntersection`

**Correction options:** `Bonferroni`, `BonferroniHolm`, `BenjaminiHochberg`, `None`

## Architecture

```mermaid
graph TD
    %% Semantic, high-contrast class definitions
    classDef file fill:#f8f9fa,stroke:#6c757d,stroke-width:2px,stroke-dasharray: 5 5,color:#212529;
    classDef struct fill:#e3f2fd,stroke:#0d6efd,stroke-width:2px,color:#052c65;
    classDef process fill:#e8f5e9,stroke:#198754,stroke-width:2px,color:#0a3622;
    classDef decision fill:#fff3cd,stroke:#ffc107,stroke-width:2px,color:#664d03;
    classDef enumType fill:#f3e5f5,stroke:#9c27b0,stroke-width:2px,color:#4a148c;

    ProblemFile[("<div style='width: 120px; text-align: center;'>Config File</div>")]:::file
    ProblemFile --> Data_Loading

    subgraph Data_Loading ["Load"]
        direction TB
  SF[("<div style='width: 120px; text-align: center;'>Study Genes</div>")]:::file
        PF[("<div style='width: 120px; text-align: center;'>Population Genes</div>")]:::file
        GOF[("<div style='width: 120px; text-align: center;'>Gene Ontology</div>")]:::file
        GAF[("<div style='width: 120px; text-align: center;'>GO Annotation</div>")]:::file

        AnnotationIndex["AnnotationIndex"]:::struct

        SF & PF & GOF & GAF --> AnnotationIndex
    end

    subgraph Analysis ["Enrichment Analysis"]
        direction TB

        subgraph Frequentist_Module ["Frequentist Approach"]
            direction TB
            subgraph FreqConfig ["Test Configuration"]
                Background["Background"]:::enumType
                Correction["Correction"]:::enumType                    
            end
            FreqTest["Test: Hypergeometric"]:::process
            FreqMeasure["Measure: P-Value"]:::struct
            
            Background & Correction --> FreqTest
            FreqTest --> FreqMeasure
        end

        subgraph Bayesian_Module ["Bayesian Approach"]
            direction TB
            subgraph BayesianData ["Data Structures"]
                Model["Model: OrModel"]:::struct
                Cache["Cache: OrCache"]:::struct
                State["State: MgsaState"]:::struct
            end

            Algorithm["Algorithm: Metropolis-Hastings"]:::process
            BayesMeasure["Measure: Posterior Probability"]:::struct

            Model & Cache & State --> Algorithm
            Algorithm --> BayesMeasure
        end
    end

    BranchDec{"Method Selection"}:::decision

    AnnotationIndex --> BranchDec

    BranchDec --> FreqConfig
    BranchDec --> BayesianData

    subgraph Output_Phase ["Output"]
        direction TB
        Result["AnalysisResult"]:::struct
        Writer["CSV Writer"]:::process
        OutputFile[("result.csv")]:::file

        FreqMeasure --> Result
        BayesMeasure --> Result

        Result --> Writer
        Writer --> OutputFile
    end
```