"""
Baseline comparison: matching criteria from simple (Gene-only) to
complete (6-field). Runs exact-set matching on the same 117-disease
pool and reports matched targets + mean candidates per target.
"""

from typing import Dict, List, Tuple
import pandas as pd

INPUT_FILE = r"C:\Users\biofa\Desktop\DATAEXTRACTION\dataset_finale_Stoichiometry_CLEAN.csv"

USECOLS = [
    "ORPHAcode", "DiseaseName", "Gene", "Entry",
    "Rhea_ID", "EC", "Cofactor", "Pathway",
    "HPO_ID", "HPO_Label",
]

FULL_PROFILE = ["Gene", "Entry", "EC", "Rhea_ID", "Cofactor", "Pathway"]

def _normalise_cofactor(val):
    if pd.isna(val):
        return val
    parts = [p.strip() for p in val.split(",") if p.strip()]
    return ", ".join(sorted(parts))

def load_and_filter() -> pd.DataFrame:
    """Load, normalise, apply cofactor/pathway filter, return clean df."""
    df = pd.read_csv(INPUT_FILE, sep=";", dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()
    for col in df.columns:
        df[col] = df[col].str.strip()
        df[col] = df[col].replace(
            {"": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NAN": pd.NA,
             "none": pd.NA, "None": pd.NA, "NONE": pd.NA,
             "null": pd.NA, "Null": pd.NA, "NULL": pd.NA}
        )
    df["Cofactor"] = df["Cofactor"].apply(_normalise_cofactor)
    # filter
    cof = set(df.loc[df["Cofactor"].notna(), "ORPHAcode"])
    pth = set(df.loc[df["Pathway"].notna(),   "ORPHAcode"])
    df  = df[df["ORPHAcode"].isin(cof & pth)].copy()
    return df


def run_matching(profile_keys: List[str], df: pd.DataFrame) -> dict:
    """
    Run exact-set matching using only `profile_keys` for profile
    construction.  The 6-field completeness requirement (dropna on
    FULL_PROFILE) is always enforced so the disease pool is identical.
    """
    # ---- HPO label map ----
    hpo_label_map: Dict[str, str] = (
        df[["HPO_ID", "HPO_Label"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates(subset=["HPO_ID"])
        .set_index("HPO_ID")["HPO_Label"]
        .to_dict()
    )

    # ---- Profiles (always complete on FULL_PROFILE) ----
    profile_df = (
        df[["ORPHAcode", "DiseaseName"] + FULL_PROFILE]
        .drop_duplicates()
        .dropna(subset=FULL_PROFILE, how="any")
    )

    # Build tuple using ONLY the requested keys
    def build_profile_set(group: pd.DataFrame) -> frozenset:
        return frozenset(
            tuple(row[k] for k in profile_keys)
            for _, row in group.iterrows()
        )

    disease_profiles = (
        profile_df[["ORPHAcode"] + profile_keys]
        .drop_duplicates()
        .groupby("ORPHAcode", group_keys=False)
        .apply(build_profile_set)
        .rename("profile_set")
        .reset_index()
    )

    disease_names = df[["ORPHAcode", "DiseaseName"]].drop_duplicates("ORPHAcode")
    disease_profiles = disease_profiles.merge(disease_names, on="ORPHAcode", how="left")

    n_profiles = len(disease_profiles)

    # ---- HPO sets ----
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
    diseases["hpo_set"] = diseases["hpo_set"].apply(
        lambda x: x if isinstance(x, set) else set()
    )
    diseases["n_hpo"] = diseases["hpo_set"].apply(len)

    # ---- Split ----
    targets = diseases[diseases["n_hpo"] == 0].reset_index(drop=True)
    donors  = diseases[diseases["n_hpo"] >  0].reset_index(drop=True)
    n_targets = len(targets)
    n_donors  = len(donors)

    # ---- Matching ----
    donor_index: Dict[frozenset, List[int]] = {}
    for idx, row in donors.iterrows():
        donor_index.setdefault(row["profile_set"], []).append(idx)

    results = []
    for _, target in targets.iterrows():
        matched = donors.loc[donor_index.get(target["profile_set"], [])]
        n_match = len(matched)
        if n_match == 0:
            results.append({"n_exact": 0, "n_candidates": 0})
            continue

        hpo_counts: Dict[str, int] = {}
        for _, donor in matched.iterrows():
            for hpo_id in donor["hpo_set"]:
                hpo_counts[hpo_id] = hpo_counts.get(hpo_id, 0) + 1

        results.append({"n_exact": n_match, "n_candidates": len(hpo_counts)})

    # ---- Stats ----
    df_res = pd.DataFrame(results)
    n_matched  = (df_res["n_exact"] > 0).sum()
    n_no_match = (df_res["n_exact"] == 0).sum()
    matched_rows = df_res[df_res["n_exact"] > 0]
    mean_candidates = matched_rows["n_candidates"].mean() if len(matched_rows) else 0
    mean_donors     = matched_rows["n_exact"].mean() if len(matched_rows) else 0

    return {
        "label":            " + ".join(profile_keys),
        "n_profiles":       n_profiles,
        "n_targets":        n_targets,
        "n_donors":         n_donors,
        "n_matched":        n_matched,
        "n_no_match":       n_no_match,
        "total_candidates": int(df_res["n_candidates"].sum()),
        "mean_candidates":  round(mean_candidates, 1),
        "mean_donors":      round(mean_donors, 1),
    }


# ===================================================================
if __name__ == "__main__":
    print("Loading and filtering data ...")
    df = load_and_filter()

    # Baseline combinations for MAIN PAPER
    MAIN_BASELINES = [
        ["Gene"],
        ["Entry"],
        ["EC"],
        ["Rhea_ID"],
        ["Pathway"],
        ["Gene", "EC", "Rhea_ID"],
        ["Gene", "Entry", "EC", "Rhea_ID"],
        ["Gene", "Entry", "EC", "Rhea_ID", "Cofactor", "Pathway"],   # 6-field
    ]

    MAIN_LABELS = [
        "Gene-only",
        "Entry-only",
        "EC-only",
        "Rhea-only",
        "Pathway-only",
        "Gene+EC+Rhea",
        "Gene+Entry+EC+Rhea",
        "Exact 6-field",
    ]

    print(f"\n{'Criterion':<30} {'Matched':>8} {'Mean cand.':>10} {'Mean donors':>11}")
    print("-" * 62)

    all_stats = []
    for keys, lbl in zip(MAIN_BASELINES, MAIN_LABELS):
        st = run_matching(keys, df)
        st["display_label"] = lbl
        all_stats.append(st)
        print(f"{lbl:<30} {st['n_matched']:>3}/{st['n_targets']}     {st['mean_candidates']:>8.1f}     {st['mean_donors']:>8.1f}")

    # ------------------------------------------------------------------
    # LaTeX table
    # ------------------------------------------------------------------
    print("\n" + "=" * 62)
    print("LATEX TABLE (main paper)")
    print("=" * 62)

    print("""
\\begin{table*}[!t]
\\centering
\\caption{Internal comparison of biologically motivated matching criteria.
         Each row reports the number of matched targets (out of 32
         unannotated targets) and the mean candidate HPO terms per
         matched target.  All criteria use the same 117-disease pool
         with complete six-field profiles.}
\\label{tab:baseline_comparison}
\\small
\\begin{tabular}{lccc}
\\toprule
\\textbf{Matching criterion} & \\textbf{Matched targets} &
\\textbf{Mean candidates} & \\textbf{Mean donors} \\\\
\\midrule""")

    for st in all_stats:
        label = st["display_label"]
        matched_str = f"{st['n_matched']}/32"
        print(f"{label:40} & {matched_str:>6}      & {st['mean_candidates']:>6.1f}        & {st['mean_donors']:>5.1f} \\\\")

    print("""\\bottomrule
\\end{tabular}
\\end{table*}""")

    # ------------------------------------------------------------------
    # Supplementary baselines
    # ------------------------------------------------------------------
    SUPP_BASELINES = [
        ["Gene", "EC"],
        ["Gene", "Rhea_ID"],
        ["EC", "Rhea_ID"],
        ["Gene", "Entry"],
        ["Gene", "Entry", "EC", "Rhea_ID", "Cofactor"],
        ["Gene", "Entry", "EC", "Rhea_ID", "Pathway"],
    ]
    SUPP_LABELS = [
        "Gene+EC",
        "Gene+Rhea",
        "EC+Rhea",
        "Gene+Entry",
        "Gene+Entry+EC+Rhea+Cofactor",
        "Gene+Entry+EC+Rhea+Pathway",
    ]

    print("\n" + "=" * 62)
    print("SUPPLEMENTARY BASELINES")
    print("=" * 62)
    print(f"{'Criterion':<40} {'Matched':>8} {'Mean cand.':>10} {'Mean donors':>11}")
    print("-" * 72)
    for keys, lbl in zip(SUPP_BASELINES, SUPP_LABELS):
        st = run_matching(keys, df)
        print(f"{lbl:<40} {st['n_matched']:>3}/{st['n_targets']}     {st['mean_candidates']:>8.1f}     {st['mean_donors']:>8.1f}")

    print("\nDone.")
