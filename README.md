# Ontologizer

Fast and safe implementation of the Ontologizer — a tool for Gene Ontology (GO) enrichment analysis using Frequentist (
hypergeometric test) and Bayesian (inference) methods.

## Gene Symbols

Gene symbols in the study and population gene sets must match the gene symbols in the `.gaf` annotation file —
specifically column 2 (`DB_Object_Symbol`).

## Quick Start

### 1. Provide a GO annotation file

Download a `.gaf` file for your organism from
the [GO Annotation Database](https://current.geneontology.org/products/pages/downloads.html) and place it in `data/`.
For example:

- Yeast: `data/goa_yeast.gaf`
- Human: `data/goa_human.gaf`

The GO ontology (`data/go-basic.json`) is downloaded automatically on first run.

### 2. Run an example

Run the yeast Bayesian example (uses `config.json` in the project root by default):

```bash
cargo run --release
```

Or point to a specific example config:

```bash
cargo run --release -- examples/yeast/config.json
cargo run --release -- examples/go0090717/config.json
```

Results are written to `output/enrichment_result.csv`.

### 3. Build the binary

```bash
cargo build --release --features cli
./target/release/ontologizer examples/yeast/config.json
```

## Configuration

Each example ships with a ready-to-use `config.json`. The schema is:

```json
{
  "study_genes_path": "examples/yeast/study_genes.txt",
  "population_genes_path": "examples/yeast/population_genes.txt",
  "go_path": "data/go-basic.json",
  "goa_path": "data/goa_yeast.gaf",
  "method": {
    "method": "bayesian"
  }
}
```

For frequentist analysis:

```json
{
  "method": {
    "method": "frequentist",
    "topology": "standard",
    "correction": "bonferroni"
  }
}
```

**Topology options:** `standard`, `parentUnion`, `parentIntersection`

**Correction options:** `bonferroni`, `bonferroniHolm`, `benjaminiHochberg`, `none`

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