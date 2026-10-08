#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage 2 exact biochemical-profile matching and HPO prioritization.

This script reproduces the final Python/Pandas analysis reported in the paper.
It excludes Maude consistency checks and starts from the integrated dataset:

    data/dataset_finale_Stoichiometry_CLEAN.csv.gz

A disease profile is the order-independent set of complete descriptors
(Gene, Entry, EC, Rhea_ID, Cofactor, Pathway). An HPO-unannotated disease is
matched to an HPO-annotated disease only when their complete profile sets are
identical.
"""

from __future__ import annotations
import argparse
from pathlib import Path
from stage2_core import (
    FULL_PROFILE,
    candidate_type_summary,
    load_integrated_dataset,
    run_exact_profile_matching,
    save_summary_table,
    target_summary,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / 'data' / 'dataset_finale_Stoichiometry_CLEAN.csv.gz'
DEFAULT_OUTPUT = REPO_ROOT / 'outputs' / 'hpo_prioritization_results.csv'

PAPER_EXPECTED_COUNTS = {
    'raw_rows': 134209,
    'raw_diseases': 1276,
    'rows_after_filter': 10609,
    'complete_profile_diseases': 117,
    'unannotated_diseases': 32,
    'annotated_diseases': 85,
    'matched_unannotated_diseases': 10,
    'unmatched_unannotated_diseases': 22,
    'candidate_associations': 495,
    'total_output_rows': 517,
    'unique_candidate_hpo_terms': 401,
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Run exact six-field biochemical-profile matching and candidate HPO prioritization.'
    )
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT,
                        help='Integrated Stage 2 CSV or CSV.GZ.')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT,
                        help='Main output CSV.')
    parser.add_argument('--check-paper-counts', action='store_true',
                        help='Fail when the input does not reproduce the paper counts.')
    return parser.parse_args()

def check_counts(dataset_stats, matching_stats) -> None:
    observed = {
        'raw_rows': dataset_stats.raw_rows,
        'raw_diseases': dataset_stats.raw_diseases,
        'rows_after_filter': dataset_stats.rows_after_filter,
        'complete_profile_diseases': matching_stats.complete_profile_diseases,
        'unannotated_diseases': matching_stats.unannotated_diseases,
        'annotated_diseases': matching_stats.annotated_diseases,
        'matched_unannotated_diseases': matching_stats.matched_unannotated_diseases,
        'unmatched_unannotated_diseases': matching_stats.unmatched_unannotated_diseases,
        'candidate_associations': matching_stats.candidate_associations,
        'total_output_rows': matching_stats.total_output_rows,
        'unique_candidate_hpo_terms': matching_stats.unique_candidate_hpo_terms,
    }
    differences = {
        name: (PAPER_EXPECTED_COUNTS[name], value)
        for name, value in observed.items()
        if value != PAPER_EXPECTED_COUNTS[name]
    }
    if differences:
        details = '\n'.join(
            f'  {name}: expected {expected}, found {found}'
            for name, (expected, found) in differences.items()
        )
        raise RuntimeError('The input does not reproduce the paper counts:\n' + details)

def main() -> None:
    args = parse_args()
    input_file = args.input if args.input.is_absolute() else REPO_ROOT / args.input
    output_file = args.output if args.output.is_absolute() else REPO_ROOT / args.output
    output_file.parent.mkdir(parents=True, exist_ok=True)

    dataframe, dataset_stats = load_integrated_dataset(
        input_file, apply_cofactor_pathway_filter=True
    )
    results, matching_stats = run_exact_profile_matching(
        dataframe,
        profile_keys=FULL_PROFILE,
        completeness_keys=FULL_PROFILE,
    )

    results.to_csv(output_file, sep=';', index=False, encoding='utf-8-sig')
    save_summary_table(dataset_stats, matching_stats,
                       output_file.parent / 'stage2_summary.csv')
    candidate_type_summary(results).to_csv(
        output_file.parent / 'candidate_type_summary.csv',
        index=False, encoding='utf-8-sig'
    )
    targets = target_summary(results)
    targets.to_csv(output_file.parent / 'target_summary.csv',
                   index=False, encoding='utf-8-sig')
    targets[targets['match_status'] == 'no_match'].to_csv(
        output_file.parent / 'unmatched_diseases.csv',
        index=False, encoding='utf-8-sig'
    )

    if args.check_paper_counts:
        check_counts(dataset_stats, matching_stats)

    print('=' * 72)
    print('STAGE 2: EXACT PROFILE MATCHING AND HPO PRIORITIZATION')
    print('=' * 72)
    print(f'Raw rows:                         {dataset_stats.raw_rows:,}')
    print(f'Raw diseases:                     {dataset_stats.raw_diseases:,}')
    print(f'Rows after cofactor/pathway filter: {dataset_stats.rows_after_filter:,}')
    print(f'Diseases with complete profiles:   {matching_stats.complete_profile_diseases:,}')
    print(f'HPO-unannotated diseases:          {matching_stats.unannotated_diseases:,}')
    print(f'HPO-annotated diseases:            {matching_stats.annotated_diseases:,}')
    print(f'Matched / unmatched:               {matching_stats.matched_unannotated_diseases}/{matching_stats.unmatched_unannotated_diseases}')
    print(f'Candidate associations:            {matching_stats.candidate_associations:,}')
    print(f'Results written to: {output_file}')

if __name__ == '__main__':
    main()
