# TODO for Carlos

Please add the technical instructions needed to run the Maude specifications.

## Needed information

1. Maude version used.
2. Tested operating system or execution environment.
3. Exact content or source of `maude/semantics/semantics.maude`.
4. Exact command for running each `*_properties.maude` file.
5. Expected output for positive reachability queries.
6. Expected output for negative reachability queries.
7. Whether any file names or module names should be adjusted before release.
8. Any known warnings, limitations, or manual steps.

## Files requiring execution instructions

- `maude/stage1_phe_tyr/phe_tyr_catabolism_properties.maude`
- `maude/stage2_consistency_checks/ampd2_purine_salvage/ampd2_purine_salvage_properties.maude`
- `maude/stage2_consistency_checks/b3galt6_glycosylation/b3galt6_glycosylation_properties.maude`
- `maude/stage2_consistency_checks/psat1_serine_biosynthesis/psat1_serine_biosynthesis_properties.maude`
- `maude/stage2_consistency_checks/sepsecs_sec_trna_biosynthesis/sepsecs_sec_trna_biosynthesis_properties.maude`

## Notes from the repository preparation

- The original uploaded Maude files used `load ./semantics.maude`.
- The split GitHub hierarchy now expects the shared semantics file under `maude/semantics/semantics.maude`.
- The raw Stage 1 Phe--Tyr source includes comments for Properties 10 and 11 but no explicit search commands after those comments. The executable P10/P11 queries were added from `paper_appendix/corrected_combined_appendices.tex`, Table XVIII.
