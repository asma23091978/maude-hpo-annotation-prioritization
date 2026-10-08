#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare exact matching criteria on the same complete-profile disease pool.

Every criterion is evaluated on diseases having complete values for all six
fields. Only the fields used to define exact profile equality are changed.
This reproduces the internal baseline comparison reported in the paper.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from stage2_core import (
    FULL_PROFILE,
    load_integrated_dataset,
    run_exact_profile_matching,
)


BASELINES = [
    ("Gene-only", ["Gene"]),
    ("UniProt Entry-only", ["Entry"]),
    ("EC-only", ["EC"]),
    ("Rhea-only", ["Rhea_ID"]),
    ("Pathway-only", ["Pathway"]),
    ("Gene+EC+Rhea", ["Gene", "EC", "Rhea_ID"]),
    (
        "Gene+Entry+EC+Rhea",
        ["Gene", "Entry", "EC", "Rhea_ID"],
    ),
    ("Exact six-field", FULL_PROFILE),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the Stage 2 internal matching baseline comparison."
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/baseline_comparison.csv"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    dataframe, _ = load_integrated_dataset(
        args.input_csv,
        apply_cofactor_pathway_filter=True,
    )

    rows = []
    for label, profile_keys in BASELINES:
        _, stats = run_exact_profile_matching(
            dataframe,
            profile_keys=profile_keys,
            # Keep the disease pool fixed across criteria.
            completeness_keys=FULL_PROFILE,
        )
        rows.append(
            {
                "criterion": label,
                "profile_fields": "+".join(profile_keys),
                "matched_diseases": stats.matched_unannotated_diseases,
                "total_unannotated_diseases": stats.unannotated_diseases,
                "mean_hpo_associations_per_matched_disease":
                    stats.mean_candidates_per_matched_disease,
                "mean_annotated_matches_per_matched_disease":
                    stats.mean_annotated_matches_per_matched_disease,
            }
        )

    output = pd.DataFrame(rows)
    output.to_csv(
        args.output,
        index=False,
        encoding="utf-8-sig",
    )

    print(output.to_string(index=False))
    print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()
