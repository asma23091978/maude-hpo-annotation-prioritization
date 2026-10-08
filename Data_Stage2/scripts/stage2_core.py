# -*- coding: utf-8 -*-
"""Core functions for Stage 2 reaction-profile matching and HPO prioritization.

This module implements the Python/Pandas analysis described in the paper's
Stage 2. It starts from the already integrated disease-reaction CSV and does
not run Maude or download data from Orphanet, UniProt, Rhea, or ChEBI.

Mathematical correspondence
---------------------------
For disease d and disease-reaction row j, a complete six-field descriptor is

    Phi_j(d) = (Gene, Entry, EC, Rhea_ID, Cofactor, Pathway).

The disease biochemical profile is the order-independent set

    P(d) = {Phi_1(d), ..., Phi_m(d)}.

The HPO annotation set is

    H(d) = {HPO_ID values associated with d}.

An HPO-unannotated disease du exactly matches an HPO-annotated disease da
when P(du) == P(da). Candidate HPO terms are collected from all exact
matches and classified by their support count.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd


REQUIRED_COLUMNS: List[str] = [
    "ORPHAcode",
    "DiseaseName",
    "Gene",
    "Entry",
    "Rhea_ID",
    "EC",
    "Cofactor",
    "Pathway",
    "HPO_ID",
    "HPO_Label",
]

FULL_PROFILE: List[str] = [
    "Gene",
    "Entry",
    "EC",
    "Rhea_ID",
    "Cofactor",
    "Pathway",
]

OUTPUT_COLUMNS: List[str] = [
    "ORPHAcode_target",
    "DiseaseName_target",
    "n_exact_donors",
    "ORPHAcode_donors",
    "DiseaseName_donors",
    "HPO_ID",
    "HPO_Label",
    "donor_count",
    "Supporting_ORPHAcode_donors",
    "candidate_type",
    "is_selected",
    "is_fully_supported",
]

MISSING_VALUE_STRINGS = {
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


@dataclass(frozen=True)
class DatasetStats:
    """Counts recorded while loading and filtering the integrated dataset."""

    raw_rows: int
    raw_diseases: int
    rows_after_filter: int
    diseases_after_filter: int


@dataclass(frozen=True)
class MatchingStats:
    """Summary statistics for one exact-profile matching run."""

    complete_profile_rows: int
    complete_profile_diseases: int
    unannotated_diseases: int
    annotated_diseases: int
    matched_unannotated_diseases: int
    unmatched_unannotated_diseases: int
    candidate_associations: int
    total_output_rows: int
    unique_candidate_hpo_terms: int
    mean_candidates_per_matched_disease: float
    median_candidates_per_matched_disease: float
    mean_annotated_matches_per_matched_disease: float

    def to_dict(self) -> dict:
        return asdict(self)


def normalise_cofactor(value: object) -> object:
    """Normalize a comma-separated cofactor list as an order-independent set.

    The original experiment treated the order of cofactors as irrelevant.
    For example, ``"Fe, BH4"`` and ``"BH4, Fe"`` both become ``"BH4, Fe"``.

    The function intentionally uses commas only, matching the final notebook
    implementation used for the paper.
    """
    if pd.isna(value):
        return pd.NA

    parts = {
        part.strip()
        for part in str(value).split(",")
        if part.strip()
    }
    return ", ".join(sorted(parts)) if parts else pd.NA


def load_integrated_dataset(
    input_file: Path | str,
    *,
    sep: str = ";",
    apply_cofactor_pathway_filter: bool = True,
) -> Tuple[pd.DataFrame, DatasetStats]:
    """Load, normalize, and optionally filter the integrated Stage 2 dataset.

    Parameters
    ----------
    input_file:
        Semicolon-separated CSV containing the columns in REQUIRED_COLUMNS.
    sep:
        CSV delimiter. The experiment used ``;``.
    apply_cofactor_pathway_filter:
        When True, retain diseases having at least one non-missing Cofactor
        value and at least one non-missing Pathway value.

    Returns
    -------
    dataframe, stats
        Cleaned row-level dataset and counts before/after the optional filter.
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input dataset not found: {input_path}")

    dataframe = pd.read_csv(
        input_path,
        sep=sep,
        dtype=str,
        low_memory=False,
    )
    dataframe.columns = dataframe.columns.str.strip()

    missing_columns = [
        column for column in REQUIRED_COLUMNS if column not in dataframe.columns
    ]
    if missing_columns:
        raise ValueError(
            "The input CSV is missing required columns: "
            + ", ".join(missing_columns)
        )

    dataframe = dataframe[REQUIRED_COLUMNS].copy()

    # Strip cell whitespace and map common textual missing values to pd.NA.
    for column in dataframe.columns:
        dataframe[column] = dataframe[column].str.strip()
        dataframe[column] = dataframe[column].replace(MISSING_VALUE_STRINGS)

    # A row without an ORPHA code cannot be assigned to a disease.
    dataframe = dataframe.dropna(subset=["ORPHAcode"]).copy()

    # Cofactor order must not affect exact profile equality.
    dataframe["Cofactor"] = dataframe["Cofactor"].apply(normalise_cofactor)

    raw_rows = len(dataframe)
    raw_diseases = dataframe["ORPHAcode"].nunique()

    if apply_cofactor_pathway_filter:
        diseases_with_cofactor = set(
            dataframe.loc[dataframe["Cofactor"].notna(), "ORPHAcode"]
        )
        diseases_with_pathway = set(
            dataframe.loc[dataframe["Pathway"].notna(), "ORPHAcode"]
        )
        retained_diseases = diseases_with_cofactor & diseases_with_pathway
        dataframe = dataframe[
            dataframe["ORPHAcode"].isin(retained_diseases)
        ].copy()

    stats = DatasetStats(
        raw_rows=raw_rows,
        raw_diseases=raw_diseases,
        rows_after_filter=len(dataframe),
        diseases_after_filter=dataframe["ORPHAcode"].nunique(),
    )
    return dataframe, stats


def _build_profile_set(group: pd.DataFrame) -> frozenset:
    """Convert one disease's reaction rows into an order-independent set."""
    return frozenset(group.itertuples(index=False, name=None))


def build_disease_table(
    dataframe: pd.DataFrame,
    *,
    profile_keys: Sequence[str] = FULL_PROFILE,
    completeness_keys: Sequence[str] = FULL_PROFILE,
) -> Tuple[pd.DataFrame, int]:
    """Build disease-level biochemical profiles P(d) and HPO sets H(d).

    ``profile_keys`` controls which fields define profile equality.
    ``completeness_keys`` controls which fields must be non-missing before a
    reaction descriptor is accepted. Baseline comparisons vary profile_keys
    but keep completeness_keys equal to all six fields, ensuring the same
    117-disease pool is used for every criterion.
    """
    unknown_profile_keys = [
        key for key in profile_keys if key not in FULL_PROFILE
    ]
    if unknown_profile_keys:
        raise ValueError(
            "Unsupported profile fields: " + ", ".join(unknown_profile_keys)
        )

    # Keep complete, non-duplicate six-field reaction descriptors.
    profile_rows = (
        dataframe[
            ["ORPHAcode", "DiseaseName", *FULL_PROFILE]
        ]
        .drop_duplicates()
        .dropna(subset=list(completeness_keys), how="any")
    )

    # Construct P(d) from only the fields requested for the current analysis.
    disease_profiles = (
        profile_rows
        .groupby("ORPHAcode", sort=False)[list(profile_keys)]
        .apply(_build_profile_set)
        .rename("profile_set")
        .reset_index()
    )

    disease_names = (
        dataframe[["ORPHAcode", "DiseaseName"]]
        .drop_duplicates(subset=["ORPHAcode"])
    )
    disease_profiles = disease_profiles.merge(
        disease_names,
        on="ORPHAcode",
        how="left",
    )

    # H(d) is built from HPO identifiers only; labels are presentation data.
    disease_hpo = (
        dataframe[["ORPHAcode", "HPO_ID"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates()
        .groupby("ORPHAcode")["HPO_ID"]
        .apply(set)
        .rename("hpo_set")
        .reset_index()
    )

    diseases = disease_profiles.merge(
        disease_hpo,
        on="ORPHAcode",
        how="left",
    )
    diseases["hpo_set"] = diseases["hpo_set"].apply(
        lambda value: value if isinstance(value, set) else set()
    )
    diseases["n_hpo"] = diseases["hpo_set"].apply(len)

    return diseases, len(profile_rows)


def run_exact_profile_matching(
    dataframe: pd.DataFrame,
    *,
    profile_keys: Sequence[str] = FULL_PROFILE,
    completeness_keys: Sequence[str] = FULL_PROFILE,
) -> Tuple[pd.DataFrame, MatchingStats]:
    """Match HPO-unannotated diseases to HPO-annotated exact profile matches.

    Candidate categories follow the final paper implementation:

    * ``single_match``: present in exactly one exact-matching annotated disease;
    * ``recurrent``: present in at least two matches, below maximum support;
    * ``dominant``: recurrent term with maximum support, but not universal;
    * ``fully_supported``: present in every exact-matching annotated disease;
    * ``no_match``: no annotated disease has the same complete profile.
    """
    hpo_label_map: Dict[str, str] = (
        dataframe[["HPO_ID", "HPO_Label"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates(subset=["HPO_ID"])
        .set_index("HPO_ID")["HPO_Label"]
        .fillna("")
        .to_dict()
    )

    diseases, complete_profile_rows = build_disease_table(
        dataframe,
        profile_keys=profile_keys,
        completeness_keys=completeness_keys,
    )

    unannotated = (
        diseases[diseases["n_hpo"] == 0]
        .reset_index(drop=True)
    )
    annotated = (
        diseases[diseases["n_hpo"] > 0]
        .reset_index(drop=True)
    )

    # Hashable frozenset profiles make exact matching an indexed lookup.
    annotated_profile_index: Dict[frozenset, List[int]] = {}
    for row_index, row in annotated.iterrows():
        annotated_profile_index.setdefault(
            row["profile_set"], []
        ).append(row_index)

    output_rows: List[dict] = []

    for _, target in unannotated.iterrows():
        matching_indices = annotated_profile_index.get(
            target["profile_set"], []
        )
        matching_annotated = annotated.loc[matching_indices]
        number_of_matches = len(matching_annotated)

        if number_of_matches == 0:
            output_rows.append(
                {
                    "ORPHAcode_target": target["ORPHAcode"],
                    "DiseaseName_target": target["DiseaseName"],
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

        # s_h(du): number of exact annotated matches carrying HPO term h.
        support_count: Dict[str, int] = {}
        supporting_orpha_codes: Dict[str, List[str]] = {}

        for _, annotated_disease in matching_annotated.iterrows():
            annotated_code = str(annotated_disease["ORPHAcode"])
            for hpo_id in annotated_disease["hpo_set"]:
                support_count[hpo_id] = support_count.get(hpo_id, 0) + 1
                supporting_orpha_codes.setdefault(hpo_id, []).append(
                    annotated_code
                )

        all_match_codes = "|".join(
            matching_annotated["ORPHAcode"].astype(str).tolist()
        )
        all_match_names = "|".join(
            matching_annotated["DiseaseName"].fillna("").astype(str).tolist()
        )

        recurrent_support = [
            count for count in support_count.values() if count >= 2
        ]
        maximum_recurrent_support = (
            max(recurrent_support) if recurrent_support else 0
        )

        for hpo_id, count in support_count.items():
            if count == 1:
                candidate_type = "single_match"
            elif count == maximum_recurrent_support:
                candidate_type = (
                    "fully_supported"
                    if count == number_of_matches
                    else "dominant"
                )
            else:
                candidate_type = "recurrent"

            # This flag is preserved from the final notebook:
            # it marks the strongest multi-match categories only.
            is_selected = candidate_type in {
                "dominant",
                "fully_supported",
            }

            output_rows.append(
                {
                    "ORPHAcode_target": target["ORPHAcode"],
                    "DiseaseName_target": target["DiseaseName"],
                    "n_exact_donors": number_of_matches,
                    "ORPHAcode_donors": all_match_codes,
                    "DiseaseName_donors": all_match_names,
                    "HPO_ID": hpo_id,
                    "HPO_Label": hpo_label_map.get(hpo_id, ""),
                    "donor_count": count,
                    "Supporting_ORPHAcode_donors": "|".join(
                        supporting_orpha_codes[hpo_id]
                    ),
                    "candidate_type": candidate_type,
                    "is_selected": is_selected,
                    "is_fully_supported": (
                        candidate_type == "fully_supported"
                    ),
                }
            )

    results = pd.DataFrame(output_rows, columns=OUTPUT_COLUMNS)

    type_order = {
        "fully_supported": 0,
        "dominant": 1,
        "recurrent": 2,
        "single_match": 3,
        "no_match": 4,
    }
    if not results.empty:
        results["_type_order"] = (
            results["candidate_type"].map(type_order).fillna(99)
        )
        results = (
            results
            .sort_values(
                [
                    "ORPHAcode_target",
                    "_type_order",
                    "donor_count",
                    "HPO_ID",
                ],
                ascending=[True, True, False, True],
                na_position="last",
            )
            .drop(columns=["_type_order"])
            .reset_index(drop=True)
        )

    candidate_rows = results[
        results["candidate_type"] != "no_match"
    ].copy()
    no_match_rows = results[
        results["candidate_type"] == "no_match"
    ].copy()

    matched_target_codes = candidate_rows["ORPHAcode_target"].drop_duplicates()
    unmatched_target_codes = no_match_rows["ORPHAcode_target"].drop_duplicates()

    candidates_per_target = (
        candidate_rows.groupby("ORPHAcode_target").size()
        if not candidate_rows.empty
        else pd.Series(dtype=float)
    )
    matches_per_target = (
        candidate_rows[
            ["ORPHAcode_target", "n_exact_donors"]
        ]
        .drop_duplicates("ORPHAcode_target")
        ["n_exact_donors"]
        if not candidate_rows.empty
        else pd.Series(dtype=float)
    )

    stats = MatchingStats(
        complete_profile_rows=complete_profile_rows,
        complete_profile_diseases=len(diseases),
        unannotated_diseases=len(unannotated),
        annotated_diseases=len(annotated),
        matched_unannotated_diseases=matched_target_codes.nunique(),
        unmatched_unannotated_diseases=unmatched_target_codes.nunique(),
        candidate_associations=len(candidate_rows),
        total_output_rows=len(results),
        unique_candidate_hpo_terms=candidate_rows["HPO_ID"].nunique(),
        mean_candidates_per_matched_disease=(
            round(float(candidates_per_target.mean()), 1)
            if not candidates_per_target.empty else 0.0
        ),
        median_candidates_per_matched_disease=(
            round(float(candidates_per_target.median()), 1)
            if not candidates_per_target.empty else 0.0
        ),
        mean_annotated_matches_per_matched_disease=(
            round(float(matches_per_target.mean()), 1)
            if not matches_per_target.empty else 0.0
        ),
    )
    return results, stats


def candidate_type_summary(results: pd.DataFrame) -> pd.DataFrame:
    """Return paper-ready counts by non-overlapping candidate category."""
    categories = [
        "single_match",
        "recurrent",
        "dominant",
        "fully_supported",
    ]
    candidate_rows = results[
        results["candidate_type"].isin(categories)
    ]

    rows = []
    for category in categories:
        subset = candidate_rows[
            candidate_rows["candidate_type"] == category
        ]
        rows.append(
            {
                "candidate_type": category,
                "target_diseases": subset["ORPHAcode_target"].nunique(),
                "unique_hpo_terms": subset["HPO_ID"].nunique(),
                "candidate_associations": len(subset),
            }
        )
    return pd.DataFrame(rows)


def target_summary(results: pd.DataFrame) -> pd.DataFrame:
    """Summarize match and candidate counts for each unannotated disease."""
    rows = []
    for (orpha_code, disease_name), group in results.groupby(
        ["ORPHAcode_target", "DiseaseName_target"],
        dropna=False,
        sort=True,
    ):
        candidate_group = group[group["candidate_type"] != "no_match"]
        rows.append(
            {
                "ORPHAcode_target": orpha_code,
                "DiseaseName_target": disease_name,
                "n_exact_donors": int(group["n_exact_donors"].max()),
                "n_candidate_associations": len(candidate_group),
                "n_single_match": int(
                    (candidate_group["candidate_type"] == "single_match").sum()
                ),
                "n_recurrent": int(
                    (candidate_group["candidate_type"] == "recurrent").sum()
                ),
                "n_dominant": int(
                    (candidate_group["candidate_type"] == "dominant").sum()
                ),
                "n_fully_supported": int(
                    (
                        candidate_group["candidate_type"]
                        == "fully_supported"
                    ).sum()
                ),
                "match_status": (
                    "matched" if len(candidate_group) else "no_match"
                ),
            }
        )
    return pd.DataFrame(rows)


def save_summary_table(
    dataset_stats: DatasetStats,
    matching_stats: MatchingStats,
    output_file: Path | str,
) -> None:
    """Write one compact two-column summary CSV."""
    rows = [
        ("raw_rows", dataset_stats.raw_rows),
        ("raw_diseases", dataset_stats.raw_diseases),
        ("rows_after_filter", dataset_stats.rows_after_filter),
        ("diseases_after_filter", dataset_stats.diseases_after_filter),
        *matching_stats.to_dict().items(),
    ]
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(
        output_file,
        index=False,
        encoding="utf-8-sig",
    )
