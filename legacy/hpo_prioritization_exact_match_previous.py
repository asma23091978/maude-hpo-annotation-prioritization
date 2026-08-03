"""
Stage 2: Rare-disease dataset expansion and HPO prioritization
==============================================================

Implements the exact reaction-profile matching method described in the paper.
The script uses repository-relative paths by default:

  input : data/dataset_finale_Stoichiometry_CLEAN.csv
  output: outputs/hpo_prioritization_results.csv

Optional command-line arguments can override these defaults:

  python src/hpo_prioritization_exact_match.py --input data/dataset_finale_Stoichiometry_CLEAN.csv --output outputs/hpo_prioritization_results.csv
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List
import argparse
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_FILE = REPO_ROOT / "data" / "dataset_finale_Stoichiometry_CLEAN.csv"
DEFAULT_OUTPUT_FILE = REPO_ROOT / "outputs" / "hpo_prioritization_results.csv"

USECOLS = [
    "ORPHAcode", "DiseaseName", "Gene", "Entry",
    "Rhea_ID", "EC", "Cofactor", "Pathway",
    "HPO_ID", "HPO_Label",
]

PROFILE_KEYS = ["Gene", "Entry", "EC", "Rhea_ID", "Cofactor", "Pathway"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Exact reaction-profile matching for candidate HPO annotation prioritization."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_FILE,
        help="Input CSV file. Default: data/dataset_finale_Stoichiometry_CLEAN.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_FILE,
        help="Output CSV file. Default: outputs/hpo_prioritization_results.csv",
    )
    return parser.parse_args()


def normalise_missing_values(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        df[col] = df[col].str.strip()
        df[col] = df[col].replace(
            {
                "": pd.NA,
                "nan": pd.NA,
                "NaN": pd.NA,
                "NAN": pd.NA,
                "none": pd.NA,
                "None": pd.NA,
                "NONE": pd.NA,
                "null": pd.NA,
                "Null": pd.NA,
                "NULL": pd.NA,
            }
        )
    return df


def normalise_cofactor(val: object) -> object:
    """Normalize cofactor strings as sorted comma-separated sets.

    This makes the comparison order-independent, so for example
    "Fe²⁺, BH4" and "BH4, Fe²⁺" are treated as the same cofactor set.
    """
    if pd.isna(val):
        return val
    parts = [p.strip() for p in str(val).split(",") if p.strip()]
    return ", ".join(sorted(parts))


def build_profile_set(group: pd.DataFrame) -> frozenset:
    """Return P(d), the set of complete six-field reaction profiles for one disease."""
    return frozenset(
        tuple(row[k] for k in PROFILE_KEYS)
        for _, row in group.iterrows()
    )


def main() -> None:
    args = parse_args()
    input_file = args.input if args.input.is_absolute() else REPO_ROOT / args.input
    output_file = args.output if args.output.is_absolute() else REPO_ROOT / args.output

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    print(f"Loading dataset from: {input_file}")
    df = pd.read_csv(input_file, sep=";", dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()

    missing_cols = sorted(set(USECOLS) - set(df.columns))
    if missing_cols:
        raise ValueError(f"Missing required columns in input file: {missing_cols}")

    df = df[USECOLS].copy()
    df = normalise_missing_values(df)
    df["Cofactor"] = df["Cofactor"].apply(normalise_cofactor)

    print(f"  Raw rows: {len(df):,}")

    # Stage 2 methodological filter:
    # retain diseases with at least one cofactor and at least one pathway annotation.
    diseases_with_cofactor = set(df.loc[df["Cofactor"].notna(), "ORPHAcode"])
    diseases_with_pathway = set(df.loc[df["Pathway"].notna(), "ORPHAcode"])
    valid_orpha_set = diseases_with_cofactor & diseases_with_pathway

    df = df[df["ORPHAcode"].isin(valid_orpha_set)].copy()
    print(f"  Rows after Stage-2 cofactor/pathway filter: {len(df):,}")

    # Canonical HPO_ID -> HPO_Label mapping: one label per HPO ID, first occurrence retained.
    hpo_label_map: Dict[str, str] = (
        df[["HPO_ID", "HPO_Label"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates(subset=["HPO_ID"])
        .set_index("HPO_ID")["HPO_Label"]
        .to_dict()
    )

    # Build reaction profiles P(d), requiring all six profile fields.
    profile_df = (
        df[["ORPHAcode", "DiseaseName"] + PROFILE_KEYS]
        .drop_duplicates()
        .dropna(subset=PROFILE_KEYS, how="any")
    )

    print("Building reaction profiles P(d) ...")
    disease_profiles = (
        profile_df
        .groupby("ORPHAcode")[PROFILE_KEYS]
        .apply(build_profile_set)
        .rename("profile_set")
        .reset_index()
    )

    disease_names = df[["ORPHAcode", "DiseaseName"]].drop_duplicates("ORPHAcode")
    disease_profiles = disease_profiles.merge(disease_names, on="ORPHAcode", how="left")

    # Build HPO sets H(d), counted strictly by HPO_ID.
    print("Building HPO sets H(d) ...")
    disease_hpo = (
        df[["ORPHAcode", "HPO_ID"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates()
        .groupby("ORPHAcode")["HPO_ID"]
        .apply(set)
        .rename("hpo_set")
        .reset_index()
    )

    diseases = disease_profiles.merge(disease_hpo, on="ORPHAcode", how="left")
    diseases["hpo_set"] = diseases["hpo_set"].apply(lambda x: x if isinstance(x, set) else set())
    diseases["n_hpo"] = diseases["hpo_set"].apply(len)

    print(f"  Total distinct valid diseases: {len(diseases):,}")

    # Classify targets and donors.
    targets = diseases[diseases["n_hpo"] == 0].reset_index(drop=True)
    donors = diseases[diseases["n_hpo"] > 0].reset_index(drop=True)

    print(f"  Target diseases  (H=∅):   {len(targets):,}")
    print(f"  Donor diseases (|H|>0):   {len(donors):,}")

    # Exact reaction-profile matching.
    print("Running exact reaction-profile matching ...")
    donor_profile_index: Dict[frozenset, List[int]] = {}
    for idx, row in donors.iterrows():
        donor_profile_index.setdefault(row["profile_set"], []).append(idx)

    results = []

    for _, target in targets.iterrows():
        d0_code = target["ORPHAcode"]
        d0_name = target["DiseaseName"]
        d0_prof = target["profile_set"]

        matched_donors = donors.loc[donor_profile_index.get(d0_prof, [])]
        n_exact = len(matched_donors)

        if n_exact == 0:
            results.append(
                {
                    "ORPHAcode_target": d0_code,
                    "DiseaseName_target": d0_name,
                    "n_exact_donors": 0,
                    "ORPHAcode_donors": "",
                    "DiseaseName_donors": "",
                    "HPO_ID": "",
                    "HPO_Label": "",
                    "donor_count": 0,
                    "Supporting_ORPHAcode_donors": "",
                    "candidate_type": "no_match",
                    "is_selected": False,
                    "is_fully_supported": False,
                }
            )
            continue

        hpo_donor_count: Dict[str, int] = {}
        hpo_donor_codes: Dict[str, List[str]] = {}

        for _, donor in matched_donors.iterrows():
            dk_code = str(donor["ORPHAcode"])
            for hpo_id in donor["hpo_set"]:
                hpo_donor_count[hpo_id] = hpo_donor_count.get(hpo_id, 0) + 1
                hpo_donor_codes.setdefault(hpo_id, []).append(dk_code)

        donor_codes_str = "|".join(matched_donors["ORPHAcode"].astype(str).tolist())
        donor_names_str = "|".join(matched_donors["DiseaseName"].fillna("").tolist())

        recurrent_counts = [c for c in hpo_donor_count.values() if c >= 2]
        max_multi_support = max(recurrent_counts) if recurrent_counts else 0

        for hpo_id, count in hpo_donor_count.items():
            if count == 1:
                ctype = "single_donor"
            elif count >= 2:
                if count == max_multi_support:
                    ctype = "fully_supported" if count == n_exact else "dominant"
                else:
                    ctype = "recurrent"
            else:
                ctype = "single_donor"

            is_selected = ctype in ("dominant", "fully_supported")
            is_fully_supported = ctype == "fully_supported"

            results.append(
                {
                    "ORPHAcode_target": d0_code,
                    "DiseaseName_target": d0_name,
                    "n_exact_donors": n_exact,
                    "ORPHAcode_donors": donor_codes_str,
                    "DiseaseName_donors": donor_names_str,
                    "HPO_ID": hpo_id,
                    "HPO_Label": hpo_label_map.get(hpo_id, ""),
                    "donor_count": count,
                    "Supporting_ORPHAcode_donors": "|".join(hpo_donor_codes[hpo_id]),
                    "candidate_type": ctype,
                    "is_selected": is_selected,
                    "is_fully_supported": is_fully_supported,
                }
            )

    # Assemble, sort, and export.
    results_df = pd.DataFrame(results)

    type_order = {
        "fully_supported": 0,
        "dominant": 1,
        "recurrent": 2,
        "single_donor": 3,
        "no_match": 4,
    }
    results_df["_type_order"] = results_df["candidate_type"].map(type_order)
    results_df = (
        results_df
        .sort_values(
            ["ORPHAcode_target", "_type_order", "donor_count", "HPO_ID"],
            ascending=[True, True, False, True],
        )
        .drop(columns=["_type_order"])
        .reset_index(drop=True)
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_file, sep=";", index=False, encoding="utf-8-sig")

    print(f"\nDone. Results written to:\n  {output_file}")
    print(f"  Total rows in output: {len(results_df):,}")

    summary = results_df.groupby("candidate_type")["ORPHAcode_target"].nunique()
    print("\nTarget diseases with >=1 candidate per type:")
    print(summary.to_string())

    matched = results_df[results_df["candidate_type"] != "no_match"]["ORPHAcode_target"].nunique()
    no_match = results_df[results_df["candidate_type"] == "no_match"]["ORPHAcode_target"].nunique()
    print(f"\nTargets with >=1 exact donor match : {matched:,}")
    print(f"Targets with no match             : {no_match:,}")


if __name__ == "__main__":
    main()
