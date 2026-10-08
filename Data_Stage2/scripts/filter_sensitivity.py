#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run the cofactor/pathway filter sensitivity analysis.

The analysis runs the same complete six-field matching pipeline twice:

* Filter ON: disease must have at least one cofactor and one pathway.
* Filter OFF: the preliminary disease-level filter is omitted.

Complete six-field reaction descriptors are still required in both runs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from stage2_core import (
    FULL_PROFILE,
    candidate_type_summary,
    load_integrated_dataset,
    run_exact_profile_matching,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Stage 2 with the cofactor/pathway filter ON/OFF."
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/filter_sensitivity"),
    )
    return parser.parse_args()


def one_run(input_csv: Path, use_filter: bool):
    dataframe, dataset_stats = load_integrated_dataset(
        input_csv,
        apply_cofactor_pathway_filter=use_filter,
    )
    results, matching_stats = run_exact_profile_matching(
        dataframe,
        profile_keys=FULL_PROFILE,
        completeness_keys=FULL_PROFILE,
    )
    return results, dataset_stats, matching_stats


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    on_results, on_data, on_match = one_run(args.input_csv, True)
    off_results, off_data, off_match = one_run(args.input_csv, False)

    on_results.to_csv(
        args.output_dir / "hpo_prioritization_filter_ON.csv",
        sep=";",
        index=False,
        encoding="utf-8-sig",
    )
    off_results.to_csv(
        args.output_dir / "hpo_prioritization_filter_OFF.csv",
        sep=";",
        index=False,
        encoding="utf-8-sig",
    )

    metrics = [
        (
            "raw_diseases",
            on_data.raw_diseases,
            off_data.raw_diseases,
        ),
        (
            "diseases_after_filter",
            on_data.diseases_after_filter,
            off_data.diseases_after_filter,
        ),
        (
            "complete_profile_diseases",
            on_match.complete_profile_diseases,
            off_match.complete_profile_diseases,
        ),
        (
            "unannotated_diseases",
            on_match.unannotated_diseases,
            off_match.unannotated_diseases,
        ),
        (
            "annotated_diseases",
            on_match.annotated_diseases,
            off_match.annotated_diseases,
        ),
        (
            "matched_unannotated_diseases",
            on_match.matched_unannotated_diseases,
            off_match.matched_unannotated_diseases,
        ),
        (
            "unmatched_unannotated_diseases",
            on_match.unmatched_unannotated_diseases,
            off_match.unmatched_unannotated_diseases,
        ),
        (
            "candidate_associations",
            on_match.candidate_associations,
            off_match.candidate_associations,
        ),
        (
            "mean_candidates_per_matched_disease",
            on_match.mean_candidates_per_matched_disease,
            off_match.mean_candidates_per_matched_disease,
        ),
        (
            "mean_annotated_matches_per_matched_disease",
            on_match.mean_annotated_matches_per_matched_disease,
            off_match.mean_annotated_matches_per_matched_disease,
        ),
    ]

    summary = pd.DataFrame(
        [
            {
                "metric": metric,
                "filter_ON": value_on,
                "filter_OFF": value_off,
                "delta_OFF_minus_ON": (
                    value_off - value_on
                    if isinstance(value_on, (int, float))
                    and isinstance(value_off, (int, float))
                    else ""
                ),
            }
            for metric, value_on, value_off in metrics
        ]
    )
    summary.to_csv(
        args.output_dir / "filter_sensitivity_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    category_on = candidate_type_summary(on_results).rename(
        columns={
            "target_diseases": "target_diseases_ON",
            "unique_hpo_terms": "unique_hpo_terms_ON",
            "candidate_associations": "candidate_associations_ON",
        }
    )
    category_off = candidate_type_summary(off_results).rename(
        columns={
            "target_diseases": "target_diseases_OFF",
            "unique_hpo_terms": "unique_hpo_terms_OFF",
            "candidate_associations": "candidate_associations_OFF",
        }
    )
    category_on.merge(
        category_off,
        on="candidate_type",
        how="outer",
    ).to_csv(
        args.output_dir / "filter_sensitivity_candidate_types.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(summary.to_string(index=False))
    print(f"\nResults written to: {args.output_dir}")


if __name__ == "__main__":
    main()
