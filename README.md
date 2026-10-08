# Reaction-Profile Matching and Maude Reachability Analysis to Prioritize Candidate Phenotype Annotations in Rare Diseases

This repository accompanies the paper *Reaction-Profile Matching and Maude Reachability Analysis to Prioritize Candidate Phenotype Annotations in Rare Diseases* (see [`paper/`](paper/)). It contains the supplied datasets, Jupyter notebooks, Maude models, property queries and generated outputs for **two scientific stages**, each organized into a data folder and a Maude folder.

## Scientific workflow

- **Stage 1 — Phe–Tyr validation:** a controlled phenylalanine–tyrosine catabolism model is used to assess qualitative biochemical reachability under normal, enzyme-perturbed and cofactor-limited conditions. The corresponding Maude property queries are supplied.
- **Stage 2 — rare-disease analysis:** biochemical reaction profiles are constructed from biological resources. Matching requires agreement on **six fields** (gene, UniProt entry, EC class, Rhea reaction, cofactor and pathway). HPO annotations are considered **after** biochemical matching to prioritize candidate terms for diseases lacking recorded HPO annotations. Separate Maude models are provided for selected biochemical consistency checks.

**Important:** Python/pandas notebooks perform the dataset preparation and exact profile matching; **Maude does not calculate the candidate HPO matches**. Maude examines reachability of biochemical entities in the encoded models. Candidate HPO terms are hypotheses for expert review, not clinically validated annotations.

## Repository structure

```text
.
├── Data_Stage1/
│   ├── Dataset_Small_validation.csv
│   └── Dataset_Small_validation.xlsx
├── Data_Stage2/
│   ├── Data_construction/
│   │   ├── Dataset_costruction.ipynb
│   │   ├── dataset_finale_Stoichiometry_CLEAN.csv
│   │   ├── dataset_finale_Stoichiometry_CLEAN.csv.gz
│   │   ├── EXTERNAL_SOURCES.md
│   │   └── ...                 # other source exports and intermediate tables
│   ├── Matching_Diseases/
│   │   ├── Diseases_matching.ipynb
│   │   ├── hpo_prioritization_results.csv
│   │   ├── stage2_supplementary_tables/
│   │   ├── reproduced_results/
│   │   └── ...                 # matching diagnostics and tables
│   └── scripts/
│       ├── Dataset_costruction.ipynb
│       └── Diseases_matching.ipynb
├── Maude_Stage1/
│   ├── syntax.maude
│   ├── semantics.maude
│   ├── phe_tyr_catabolism_model.maude
│   ├── phe_tyr_catabolism_properties.maude
│   └── original_sources/
├── Maude_Stage2/
│   ├── syntax.maude
│   ├── semantics.maude
│   ├── ampd2_purine_salvage/
│   ├── b3galt6_glycosylation/
│   ├── psat1_serine_biosynthesis/
│   ├── sepsecs_sec_trna_biosynthesis/
│   └── original_sources/
├── paper/
├── requirements.txt
└── README.md
```

The files in `Data_Stage2/scripts/` are **copies of the two notebooks**, not standalone `.py` programs. The folders `original_sources/` retain the supplied original Maude text files alongside the organized `.maude` models.

## Requirements

### Maude

The intended version is **Maude 3.5**, or a compatible version from the 3.5 series. Maude is installed **separately** from Python; it cannot be installed by adding a GitHub link to `requirements.txt`.

- Official project: https://github.com/maude-lang/Maude
- Releases and installation information: https://github.com/maude-lang/Maude/releases

The Maude files also contain `load model-checker .`, which relies on the corresponding Maude model-checking module being available in the installation.

### Python and Jupyter

For the Python data analysis, use Python 3 with `pandas` as specified in [`requirements.txt`](requirements.txt):

```bash
python -m pip install -r requirements.txt
```

For **running the supplied Jupyter notebooks**, additional libraries used in their code are `numpy`, `matplotlib`, `lxml`, `requests` and `seaborn`, as well as Jupyter itself. These notebook-specific dependencies are **not all recorded** in `requirements.txt`; you can install them separately, for example:

```bash
python -m pip install notebook numpy matplotlib lxml requests seaborn
```

The notebooks have not been verified to run from start to finish in a clean environment. Some cells download data from external resources and depend on local file paths.

## Running Stage 1 in Maude

**Run the following commands from the repository root** (the directory containing this README). The `.maude` source files use repository-root-relative `load` paths: for example, `load Maude_Stage1/semantics.maude`, and the latter loads `Maude_Stage1/syntax.maude`.

The existing `phe_tyr_catabolism_properties.maude` file contains the reachability searches; **it is the reference for the Stage 1 execution commands**:

```bash
maude Maude_Stage1/phe_tyr_catabolism_properties.maude
```

The underlying model is in `Maude_Stage1/phe_tyr_catabolism_model.maude`. Do not invoke the properties file from inside `Maude_Stage1/` without adapting the relative `load` paths.

## Running Stage 2 Maude consistency checks

Four selected reaction-level models and their own property-query files are included. From the repository root, run each property file separately:

```bash
maude Maude_Stage2/ampd2_purine_salvage/ampd2_purine_salvage_properties.maude
maude Maude_Stage2/b3galt6_glycosylation/b3galt6_glycosylation_properties.maude
maude Maude_Stage2/psat1_serine_biosynthesis/psat1_serine_biosynthesis_properties.maude
maude Maude_Stage2/sepsecs_sec_trna_biosynthesis/sepsecs_sec_trna_biosynthesis_properties.maude
```

Each property file loads its corresponding model; the model loads `Maude_Stage2/semantics.maude`, which loads `Maude_Stage2/syntax.maude`.

**Validation status:** the presence of these scripts and valid local file paths does not establish that all Maude searches pass. Their outcomes should be checked with an installed Maude interpreter before making a full reproducibility claim. A completed negative search (`No solution`) should be distinguished from an unfinished or interrupted search.

## Stage 2 data and matching notebooks

- `Data_Stage2/Data_construction/Dataset_costruction.ipynb` contains the data-construction workflow, with raw/source data and intermediate CSV/TSV/XML files in the same `Data_construction/` folder.
- `Data_Stage2/Matching_Diseases/Diseases_matching.ipynb` contains matching and candidate-prioritization analyses; the folder also contains saved output tables and supplementary material.
- The consolidated dataset is `Data_Stage2/Data_construction/dataset_finale_Stoichiometry_CLEAN.csv`; a compressed copy with suffix `.csv.gz` is also available **in that same folder**.

**Notebook path caveat:** some notebook cells use relative paths such as `data/dataset_finale_Stoichiometry_CLEAN.csv` and `data/en_product6.xml`. Such paths do **not** match the four-folder layout automatically. Before running those cells, either update their input paths to `Data_Stage2/Data_construction/` or prepare the expected `data/` working folder. Do not interpret the saved notebook outputs as a fresh, verified end-to-end run.

The paper reports **1,276 diseases**, including **117 diseases with complete six-field profiles** (32 HPO-unannotated and 85 annotated). Exact matching identified annotated matches for **10 of the 32 HPO-unannotated diseases**, yielding **495 candidate disease–HPO associations**. These are the values to compare against the saved matching results, using the paper's described filtering and aggregation rules.

