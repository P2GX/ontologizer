# ontologizer
Fast and safe implementation of the Ontologizer 



To run the example program, enter
```bash
 cargo run --bin onto --features="cli" 
```
For faster performance, enter
```bash
cargo run --release --bin onto --features="cli"
```



## To build the binary demo (with clap)
```bash
cargo build --release --features cli
```
(the binary is then in ``./target/release/rpt``)
to run it
```bash
cargo run --features cli --bin rpt
```
## To see private features in documentation
```bash
cargo doc --document-private-items --open
```
## Structure
```mermaid
graph TD;
    %% --- Define specific styles for clarity ---
    classDef file fill:#e1ecf4,stroke:#7a99ac,stroke-width:2px,stroke-dasharray: 5 5;
    classDef struct fill:#f9f2f4,stroke:#c07b8f,stroke-width:2px;
    classDef process fill:#fff7d1,stroke:#dccc73,stroke-width:2px;

    %% --- Main Entry ---
    %% ProblemFile now acts as the controller for both Data Loading and the Method Decision
    ProblemFile[("Problem File")]:::file 
    ProblemFile --> Data_Loading
    
    %% --- FIX 1: Explicitly show Problem File controls the Method ---
    %% Dashed line implies "Configuration/Control" rather than data transformation
    ProblemFile -.->|Specifies Method| BranchDec{"Method?"}

    %% --- Data Loading Phase ---
    subgraph Data_Loading ["Data Loading"]
        direction TB
        SF[("Study Genes")]:::file
        PF[("Population Genes")]:::file
        GOF[("Gene Ontology")]:::file
        GAF[("GO Annotation")]:::file

        AnnotationIndex["Annotation Index"]:::struct
        AnnotationMap["Annotation Map"]:::struct
        
        SF & PF & GOF & GAF --> AnnotationIndex
        SF & PF & GOF & GAF --> AnnotationMap
    end

    %% --- Enrichment Analysis Block ---
    subgraph Analysis ["Enrichment Analysis"]
        direction TB
        
        %% --- BRANCH 1: Frequentist ---
        subgraph Frequentist_Module ["Frequentist"]
            direction TB
            FreqData["Data Structure"]:::struct
            FreqTest["Test:<br/>Hypergeometric"]:::process
            FreqMeasure["Measure: P-Value<br/>(Floats)"]:::struct
            
            FreqData --> FreqTest
            FreqTest -->|Produces| FreqMeasure
        end

        %% --- BRANCH 2: Bayesian ---
        subgraph Bayesian_Module ["Bayesian"]
            direction TB
            subgraph BayesianData["Data Structure"]
                Model["Model: OrModel"]:::struct
                Cache["Cache: OrCache"]:::struct
                State["State: State"]:::struct
            end

            Algorithm["Algorithm:<br/>Metropolis-Hasting"]:::process
            BayesMeasure["Measure: Probability<br/>(Floats)"]:::struct

            %% Algorithm and State create the measure
            Model & Cache & State --> Algorithm
            Algorithm -->|Create| BayesMeasure
        end
    end

    %% --- FIX 2: Data Flow (Derivation) ---
    %% Show that the Analysis Structs are DERIVED FROM the Index/Map
    AnnotationIndex & AnnotationMap --> FreqData
    AnnotationIndex & AnnotationMap --> BayesianData

    %% --- Control Flow (Decision) ---
    %% The decision activates the module (Control), separate from the data feeding it
    BranchDec -.-> Frequentist_Module
    BranchDec -.-> Bayesian_Module

    %% --- Convergence and Output ---
    subgraph Output_Phase ["Output Generation"]
        Result["EnrichmentResult"]:::process
        Writer["CSV Writer"]:::process
        OutputFile[("Final Output.csv")]:::file

        %% Both branches feed into normalization
        FreqMeasure --> Result
        BayesMeasure --> Result
        
        Result --> Writer
        Writer --> OutputFile
    end
```
