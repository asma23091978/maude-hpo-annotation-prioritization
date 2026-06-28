# Maude models and reachability properties

The original uploaded Maude `.txt` files have been preserved in `original_sources/`.

For GitHub readability and reproducibility, each source was split into two files:

- a `*_model.maude` file containing the model, entities, reactions, environments, and input processes;
- a `*_properties.maude` file containing the reachability/model-checking queries.

The shared `semantics.maude` file is still required and must be added under `maude/semantics/`.

## Split summary

| Model | Search queries | Model-checking queries |
|---|---:|---:|
| `phe_tyr_catabolism` | 47 | 3 |
| `ampd2_purine_salvage` | 8 | 0 |
| `b3galt6_glycosylation` | 7 | 0 |
| `psat1_serine_biosynthesis` | 9 | 0 |
| `sepsecs_sec_trna_biosynthesis` | 8 | 0 |


Note: the raw Stage 1 Phe--Tyr source contained descriptive comments for Properties 10 and 11 but no explicit search commands after those comments. The missing P10/P11 executable queries were added from `paper_appendix/corrected_combined_appendices.tex`, Table XVIII, not invented.
