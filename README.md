# Reaction-profile matching and Maude analysis for rare enzymatic diseases

This repository contains the data, Python code, supplementary outputs, and
Maude models associated with the reaction-centered workflow for prioritizing
candidate Human Phenotype Ontology (HPO) annotations in rare enzymatic
diseases.

The workflow has two parts:

1. **Stage 1** validates a qualitative phenylalanine-tyrosine pathway model
   through Maude reachability analysis.
2. **Stage 2** builds complete biochemical reaction profiles and matches
   HPO-unannotated diseases to HPO-annotated diseases by exact identity of six
   fields: causal gene, reviewed UniProt entry, EC class, Rhea reaction,
   cofactor, and pathway.

Maude is not used to calculate the large-scale exact profile matches. The
matching and HPO prioritization analysis is implemented in Python with Pandas.

## Repository structure

```text
.
├── data/
│   └── dataset_finale_Stoichiometry_CLEAN.csv.gz
├── src/
│   ├── stage2_core.py
│   ├── hpo_prioritization_exact_match.py
│   ├── filter_sensitivity.py
│   └── baseline_comparison.py
├── notebooks/
│   ├── enzyme_desease_extract.ipynb
│   └── final_paper_cells/
├── outputs/
├── stage1_phe_tyr/
├── stage2_consistency_checks/
├── semantics/
├── original_sources/
├── paper_appendix/
├── notes/
└── legacy/
```

## Python installation

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
pip install -r requirements.txt
```

Linux or macOS:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the final Stage 2 analysis

From the repository root:

```bash
python src/hpo_prioritization_exact_match.py --check-paper-counts
```

The main output is `outputs/hpo_prioritization_results.csv`. The script also
generates summary, category, target, and unmatched-disease tables in
`outputs/`.

## Filter-sensitivity analysis

```bash
python src/filter_sensitivity.py \
    data/dataset_finale_Stoichiometry_CLEAN.csv.gz \
    --output-dir outputs/filter_sensitivity
```

## Internal baseline comparison

```bash
python src/baseline_comparison.py \
    data/dataset_finale_Stoichiometry_CLEAN.csv.gz \
    --output outputs/baseline_comparison.csv
```

## Expected Stage 2 results

| Metric | Value |
|---|---:|
| Disease-reaction rows | 134,209 |
| Rare diseases | 1,276 |
| Rows after cofactor/pathway filtering | 10,609 |
| Diseases with complete six-field profiles | 117 |
| HPO-unannotated diseases | 32 |
| HPO-annotated diseases | 85 |
| Diseases with an exact match | 10 |
| Diseases without an exact match | 22 |
| Candidate disease-HPO associations | 495 |
| Unique candidate HPO terms | 401 |

Candidate categories:

| Category | Associations |
|---|---:|
| Single-match | 478 |
| Recurrent | 16 |
| Dominant | 1 |
| Fully supported | 0 |

## Relation to the original notebook

The final paper results are traced to cells 152, 154, and 156. Exact copies
are stored in `notebooks/final_paper_cells/`. The cleaned implementations are
in `src/`. The previous script is preserved in `legacy/`; it uses the older
term `single_donor`, whereas the final paper uses `single_match`.

## Maude dependency

The shared file `semantics/semantics.maude` is not present in the downloaded
repository snapshot. Running the Maude models requires the exact shared
semantics file used in the experiments.
