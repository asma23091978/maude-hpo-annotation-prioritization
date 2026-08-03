# -*- coding: utf-8 -*-
"""
Stage 2: Rare-disease dataset expansion and HPO prioritization
==============================================================

Implements exact reaction-profile matching between:

    unannotated_diseases:
        diseases with a complete biochemical profile but no HPO annotation;

    annotated_diseases:
        diseases with the same type of biochemical profile and at least
        one HPO annotation.

For every unannotated disease, the script identifies annotated diseases
having an identical reaction profile and prioritizes their HPO terms.

Compatibility
-------------
This version does not use:

    include_groups=False

and is therefore compatible with pandas versions older than 2.2.

Corrections applied
-------------------
A. CSV column names are stripped before selecting USECOLS.
B. Every biochemical profile tuple must contain all six fields.
C. HPO terms are counted using HPO_ID only.
D. Supporting ORPHA codes are exported for each HPO candidate.
E. Selection and full-support flags are exported.
F. Cofactors are normalized as order-independent sets.
G. Naming:
       targets       -> unannotated_diseases
       donors        -> annotated_diseases
       single_donor  -> single_match
"""

from typing import Dict, List

import pandas as pd


# ===========================================================================
# 0. CONFIGURATION
# ===========================================================================

INPUT_FILE = (
    r"C:\Users\biofa\Desktop\DATAEXTRACTION"
    r"\dataset_finale_Stoichiometry_CLEAN.csv"
)

OUTPUT_FILE = (
    r"C:\Users\biofa\Desktop\DATAEXTRACTION"
    r"\hpo_prioritization_results.csv"
)


USECOLS = [
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


# Ordered fields defining one biochemical reaction tuple Phi_j(d)
PROFILE_KEYS = [
    "Gene",
    "Entry",
    "EC",
    "Rhea_ID",
    "Cofactor",
    "Pathway",
]


# Columns exported in the final CSV.
# Declaring them explicitly also prevents errors if no result is generated.
OUTPUT_COLUMNS = [
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


# ===========================================================================
# 1. HELPER FUNCTIONS
# ===========================================================================

def normalise_cofactor(value):
    """
    Normalize a cofactor field as an order-independent set.

    Examples
    --------
    "Fe²⁺, BH4" -> "BH4, Fe²⁺"
    "BH4, Fe²⁺" -> "BH4, Fe²⁺"
    """
    if pd.isna(value):
        return value

    parts = [
        part.strip()
        for part in str(value).split(",")
        if part.strip()
    ]

    # set() also removes accidental repeated cofactors.
    return ", ".join(sorted(set(parts)))


def build_profile_set(group: pd.DataFrame) -> frozenset:
    """
    Build the biochemical profile P(d) for one disease.

    The input contains only PROFILE_KEYS because the columns are explicitly
    selected before GroupBy.apply(). Therefore, ORPHAcode is never included
    in the biochemical tuple.

    Returns
    -------
    frozenset
        Set of complete six-field biochemical reaction tuples.
    """
    return frozenset(
        group.itertuples(index=False, name=None)
    )


# ===========================================================================
# 2. LOAD AND CLEAN THE DATASET
# ===========================================================================

print("=" * 72)
print("STAGE 2: REACTION-PROFILE MATCHING AND HPO PRIORITIZATION")
print("=" * 72)

print("\nLoading dataset ...")

df = pd.read_csv(
    INPUT_FILE,
    sep=";",
    dtype=str,
    low_memory=False,
)


# Strip spaces from column names before selecting USECOLS.
df.columns = df.columns.str.strip()


# Check that all required columns exist.
missing_columns = [
    column
    for column in USECOLS
    if column not in df.columns
]

if missing_columns:
    raise ValueError(
        "The input file is missing the following required columns: "
        + ", ".join(missing_columns)
    )


df = df[USECOLS].copy()


# Normalize whitespace and missing-value representations.
missing_value_strings = {
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


for column in df.columns:
    df[column] = df[column].str.strip()
    df[column] = df[column].replace(missing_value_strings)


# Rows without an ORPHA code cannot be assigned to a disease.
df = df.dropna(subset=["ORPHAcode"]).copy()


# Normalize cofactors before profile construction.
df["Cofactor"] = df["Cofactor"].apply(normalise_cofactor)


print("  Raw valid rows: %d" % len(df))
print("  Raw distinct diseases: %d" % df["ORPHAcode"].nunique())


# ===========================================================================
# 3. STAGE-2 COFACTOR/PATHWAY FILTER
# ===========================================================================

print("\nApplying the Stage-2 cofactor/pathway filter ...")


# A disease is retained only if it has:
#   1. at least one non-null Cofactor value;
#   2. at least one non-null Pathway value.
diseases_with_cofactor = set(
    df.loc[
        df["Cofactor"].notna(),
        "ORPHAcode",
    ]
)

diseases_with_pathway = set(
    df.loc[
        df["Pathway"].notna(),
        "ORPHAcode",
    ]
)

valid_orpha_set = (
    diseases_with_cofactor
    & diseases_with_pathway
)


df = df[
    df["ORPHAcode"].isin(valid_orpha_set)
].copy()


print(
    "  Rows after cofactor/pathway filter: %d"
    % len(df)
)

print(
    "  Diseases after cofactor/pathway filter: %d"
    % df["ORPHAcode"].nunique()
)


# ===========================================================================
# 4. BUILD THE CANONICAL HPO LABEL MAP
# ===========================================================================

# One canonical label is retained for every HPO identifier.
# If textual label variants exist, the first available value is used.
hpo_label_map: Dict[str, str] = (
    df[
        ["HPO_ID", "HPO_Label"]
    ]
    .dropna(subset=["HPO_ID"])
    .drop_duplicates(subset=["HPO_ID"])
    .set_index("HPO_ID")["HPO_Label"]
    .fillna("")
    .to_dict()
)


# ===========================================================================
# 5. BUILD BIOCHEMICAL REACTION PROFILES P(d)
# ===========================================================================

print("\nBuilding biochemical reaction profiles P(d) ...")


# Keep only complete six-field reaction tuples.
profile_df = (
    df[
        ["ORPHAcode", "DiseaseName"]
        + PROFILE_KEYS
    ]
    .drop_duplicates()
    .dropna(
        subset=PROFILE_KEYS,
        how="any",
    )
)


# ---------------------------------------------------------------------------
# IMPORTANT PANDAS COMPATIBILITY FIX
#
# PROFILE_KEYS are selected before apply().
# Therefore build_profile_set() receives only the six biochemical columns.
#
# No include_groups=False argument is used.
# ---------------------------------------------------------------------------

disease_profiles = (
    profile_df
    .groupby(
        "ORPHAcode",
        sort=False,
    )[PROFILE_KEYS]
    .apply(build_profile_set)
    .rename("profile_set")
    .reset_index()
)


# Add one disease name for each ORPHA code.
disease_names = (
    df[
        ["ORPHAcode", "DiseaseName"]
    ]
    .drop_duplicates(subset=["ORPHAcode"])
)


disease_profiles = disease_profiles.merge(
    disease_names,
    on="ORPHAcode",
    how="left",
)


print(
    "  Diseases with a complete biochemical profile: %d"
    % len(disease_profiles)
)


# ===========================================================================
# 6. BUILD HPO SETS H(d)
# ===========================================================================

print("\nBuilding HPO annotation sets H(d) ...")


# HPO terms are represented by identifiers only.
# Labels are added later using hpo_label_map.
disease_hpo = (
    df[
        ["ORPHAcode", "HPO_ID"]
    ]
    .dropna(subset=["HPO_ID"])
    .drop_duplicates()
    .groupby("ORPHAcode")["HPO_ID"]
    .apply(set)
    .rename("hpo_set")
    .reset_index()
)


# Keep every disease that has a complete biochemical profile,
# including diseases without HPO annotations.
diseases = disease_profiles.merge(
    disease_hpo,
    on="ORPHAcode",
    how="left",
)


# Diseases without HPO annotations receive an empty set.
diseases["hpo_set"] = diseases["hpo_set"].apply(
    lambda value: value
    if isinstance(value, set)
    else set()
)


diseases["n_hpo"] = diseases["hpo_set"].apply(len)


print(
    "  Total valid diseases: %d"
    % len(diseases)
)


# ===========================================================================
# 7. CLASSIFY UNANNOTATED AND ANNOTATED DISEASES
# ===========================================================================

# Unannotated diseases:
# complete biochemical profile, but H(d) is empty.
unannotated_diseases = (
    diseases[
        diseases["n_hpo"] == 0
    ]
    .reset_index(drop=True)
)


# Annotated diseases:
# complete biochemical profile and at least one HPO term.
annotated_diseases = (
    diseases[
        diseases["n_hpo"] > 0
    ]
    .reset_index(drop=True)
)


print(
    "  Unannotated diseases, H(d)=empty: %d"
    % len(unannotated_diseases)
)

print(
    "  Annotated diseases, |H(d)|>0:   %d"
    % len(annotated_diseases)
)


# ===========================================================================
# 8. INDEX ANNOTATED DISEASES BY EXACT BIOCHEMICAL PROFILE
# ===========================================================================

print("\nIndexing annotated diseases by biochemical profile ...")


# Dictionary structure:
#
#     biochemical profile -> row indices in annotated_diseases
#
# Because profile_set is a frozenset, it can be used as a dictionary key.
annotated_profile_index: Dict[frozenset, List[int]] = {}


for row_index, row in annotated_diseases.iterrows():
    profile = row["profile_set"]

    annotated_profile_index.setdefault(
        profile,
        [],
    ).append(row_index)


print(
    "  Distinct annotated biochemical profiles: %d"
    % len(annotated_profile_index)
)


# ===========================================================================
# 9. EXACT PROFILE MATCHING AND HPO PRIORITIZATION
# ===========================================================================

print("\nRunning exact biochemical-profile matching ...")


results = []


for _, unannotated_disease in unannotated_diseases.iterrows():

    target_orpha_code = unannotated_disease["ORPHAcode"]
    target_name = unannotated_disease["DiseaseName"]
    target_profile = unannotated_disease["profile_set"]


    # Retrieve annotated diseases having exactly the same frozenset profile.
    matching_indices = annotated_profile_index.get(
        target_profile,
        [],
    )


    matched_annotated_diseases = annotated_diseases.loc[
        matching_indices
    ]


    n_exact_matches = len(
        matched_annotated_diseases
    )


    # -----------------------------------------------------------------------
    # No annotated disease has the same biochemical profile.
    # -----------------------------------------------------------------------

    if n_exact_matches == 0:

        results.append(
            {
                "ORPHAcode_target": target_orpha_code,
                "DiseaseName_target": target_name,
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


    # -----------------------------------------------------------------------
    # Count support for each HPO term.
    #
    # hpo_match_count:
    #     number of matched annotated diseases carrying the HPO term.
    #
    # hpo_matching_orpha_codes:
    #     exact annotated diseases supporting that HPO term.
    # -----------------------------------------------------------------------

    hpo_match_count: Dict[str, int] = {}

    hpo_matching_orpha_codes: Dict[
        str,
        List[str],
    ] = {}


    for _, annotated_disease in matched_annotated_diseases.iterrows():

        annotated_orpha_code = str(
            annotated_disease["ORPHAcode"]
        )


        for hpo_id in annotated_disease["hpo_set"]:

            hpo_match_count[hpo_id] = (
                hpo_match_count.get(hpo_id, 0)
                + 1
            )

            hpo_matching_orpha_codes.setdefault(
                hpo_id,
                [],
            ).append(annotated_orpha_code)


    # All exact matched annotated diseases for this target.
    matched_orpha_codes_string = "|".join(
        matched_annotated_diseases[
            "ORPHAcode"
        ]
        .astype(str)
        .tolist()
    )


    matched_disease_names_string = "|".join(
        matched_annotated_diseases[
            "DiseaseName"
        ]
        .fillna("")
        .astype(str)
        .tolist()
    )


    # Maximum support among terms occurring in at least two matches.
    recurrent_support_values = [
        count
        for count in hpo_match_count.values()
        if count >= 2
    ]


    maximum_recurrent_support = (
        max(recurrent_support_values)
        if recurrent_support_values
        else 0
    )


    # -----------------------------------------------------------------------
    # Candidate classification
    # -----------------------------------------------------------------------

    for hpo_id, support_count in hpo_match_count.items():

        if support_count == 1:

            candidate_type = "single_match"


        elif support_count >= 2:

            if support_count == maximum_recurrent_support:

                if support_count == n_exact_matches:
                    candidate_type = "fully_supported"

                else:
                    candidate_type = "dominant"

            else:
                candidate_type = "recurrent"


        else:
            # Defensive fallback; normally impossible.
            candidate_type = "single_match"


        is_selected = candidate_type in {
            "dominant",
            "fully_supported",
        }


        is_fully_supported = (
            candidate_type == "fully_supported"
        )


        results.append(
            {
                "ORPHAcode_target":
                    target_orpha_code,

                "DiseaseName_target":
                    target_name,

                "n_exact_donors":
                    n_exact_matches,

                "ORPHAcode_donors":
                    matched_orpha_codes_string,

                "DiseaseName_donors":
                    matched_disease_names_string,

                "HPO_ID":
                    hpo_id,

                "HPO_Label":
                    hpo_label_map.get(hpo_id, ""),

                "donor_count":
                    support_count,

                "Supporting_ORPHAcode_donors":
                    "|".join(
                        hpo_matching_orpha_codes[hpo_id]
                    ),

                "candidate_type":
                    candidate_type,

                "is_selected":
                    is_selected,

                "is_fully_supported":
                    is_fully_supported,
            }
        )


# ===========================================================================
# 10. ASSEMBLE AND SORT THE OUTPUT
# ===========================================================================

results_df = pd.DataFrame(
    results,
    columns=OUTPUT_COLUMNS,
)


type_order = {
    "fully_supported": 0,
    "dominant": 1,
    "recurrent": 2,
    "single_match": 3,
    "no_match": 4,
}


if not results_df.empty:

    results_df["_type_order"] = (
        results_df["candidate_type"]
        .map(type_order)
        .fillna(99)
    )


    results_df = (
        results_df
        .sort_values(
            by=[
                "ORPHAcode_target",
                "_type_order",
                "donor_count",
                "HPO_ID",
            ],
            ascending=[
                True,
                True,
                False,
                True,
            ],
            na_position="last",
        )
        .drop(columns=["_type_order"])
        .reset_index(drop=True)
    )


# ===========================================================================
# 11. EXPORT RESULTS
# ===========================================================================

results_df.to_csv(
    OUTPUT_FILE,
    sep=";",
    index=False,
    encoding="utf-8-sig",
)


print("\nResults written to:")
print("  %s" % OUTPUT_FILE)

print(
    "  Total output rows: %d"
    % len(results_df)
)


# ===========================================================================
# 12. SUMMARY
# ===========================================================================

print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)


print(
    "Unannotated diseases examined: %d"
    % len(unannotated_diseases)
)

print(
    "Annotated diseases available:  %d"
    % len(annotated_diseases)
)


if results_df.empty:

    print(
        "\nNo output rows were generated. "
        "Check whether the dataset contains diseases "
        "with complete biochemical profiles."
    )

else:

    summary = (
        results_df
        .groupby("candidate_type")[
            "ORPHAcode_target"
        ]
        .nunique()
        .reindex(
            [
                "fully_supported",
                "dominant",
                "recurrent",
                "single_match",
                "no_match",
            ],
            fill_value=0,
        )
    )


    print(
        "\nNumber of unannotated diseases "
        "with candidates in each class:"
    )

    print(summary.to_string())


    matched_targets = (
        results_df.loc[
            results_df["candidate_type"] != "no_match",
            "ORPHAcode_target",
        ]
        .nunique()
    )


    unmatched_targets = (
        results_df.loc[
            results_df["candidate_type"] == "no_match",
            "ORPHAcode_target",
        ]
        .nunique()
    )


    selected_candidates = (
        results_df[
            results_df["is_selected"] == True
        ]
        .shape[0]
    )


    print(
        "\nUnannotated diseases with >=1 exact match: %d"
        % matched_targets
    )

    print(
        "Unannotated diseases with no exact match:  %d"
        % unmatched_targets
    )

    print(
        "Selected dominant/fully-supported HPO rows: %d"
        % selected_candidates
    )


print("\nDone.")