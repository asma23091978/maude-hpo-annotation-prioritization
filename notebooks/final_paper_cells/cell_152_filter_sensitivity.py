"""
Sensitivity Analysis: cofactor/pathway filter ON vs OFF
=======================================================
Compares the Stage-2 pipeline with and without the methodological
filter that requires at least one Cofactor AND at least one Pathway
per disease.
"""

from typing import Dict, List, Tuple
import pandas as pd
import sys

INPUT_FILE = r"C:\Users\biofa\Desktop\DATAEXTRACTION\dataset_finale_Stoichiometry_CLEAN.csv"

USECOLS = [
    "ORPHAcode", "DiseaseName", "Gene", "Entry",
    "Rhea_ID", "EC", "Cofactor", "Pathway",
    "HPO_ID", "HPO_Label",
]

PROFILE_KEYS = ["Gene", "Entry", "EC", "Rhea_ID", "Cofactor", "Pathway"]


def _normalise_cofactor(val):
    if pd.isna(val):
        return val
    parts = [p.strip() for p in val.split(",") if p.strip()]
    return ", ".join(sorted(parts))


def run_pipeline(use_cofactor_pathway_filter: bool) -> Tuple[pd.DataFrame, dict]:
    """
    Run the Stage-2 HPO prioritisation pipeline.

    Parameters
    ----------
    use_cofactor_pathway_filter : bool
        If True, only diseases with >=1 Cofactor AND >=1 Pathway are kept.
        If False, all diseases are retained (no methodological filter).

    Returns
    -------
    results_df : pd.DataFrame
        Full prioritisation output (same schema as hpo_prioritization_results.csv).
    stats : dict
        Summary statistics for comparison.
    """
    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    df = pd.read_csv(INPUT_FILE, sep=";", dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()
    df = df[USECOLS].copy()

    for col in df.columns:
        df[col] = df[col].str.strip()
        df[col] = df[col].replace(
            {
                "": pd.NA, "nan": pd.NA, "NaN": pd.NA, "NAN": pd.NA,
                "none": pd.NA, "None": pd.NA, "NONE": pd.NA,
                "null": pd.NA, "Null": pd.NA, "NULL": pd.NA,
            }
        )

    df["Cofactor"] = df["Cofactor"].apply(_normalise_cofactor)

    raw_rows = len(df)
    raw_diseases = df["ORPHAcode"].nunique()
    raw_diseases_with_hpo = df.dropna(subset=["HPO_ID"])["ORPHAcode"].nunique()

    # ------------------------------------------------------------------
    # Stage-2 methodological filter (optional)
    # ------------------------------------------------------------------
    if use_cofactor_pathway_filter:
        diseases_with_cofactor = set(df.loc[df["Cofactor"].notna(), "ORPHAcode"])
        diseases_with_pathway  = set(df.loc[df["Pathway"].notna(), "ORPHAcode"])
        valid_orpha_set        = diseases_with_cofactor & diseases_with_pathway
        df = df[df["ORPHAcode"].isin(valid_orpha_set)].copy()

    rows_after_filter = len(df)
    diseases_after_filter = df["ORPHAcode"].nunique()

    # ------------------------------------------------------------------
    # HPO label map
    # ------------------------------------------------------------------
    hpo_label_map: Dict[str, str] = (
        df[["HPO_ID", "HPO_Label"]]
        .dropna(subset=["HPO_ID"])
        .drop_duplicates(subset=["HPO_ID"])
        .set_index("HPO_ID")["HPO_Label"]
        .to_dict()
    )

    # ------------------------------------------------------------------
    # 2. Build reaction profiles P(d)  — requires all 6 fields
    # ------------------------------------------------------------------
    profile_df = (
        df[["ORPHAcode", "DiseaseName"] + PROFILE_KEYS]
        .drop_duplicates()
        .dropna(subset=PROFILE_KEYS, how="any")
    )

    def build_profile_set(group: pd.DataFrame) -> frozenset:
        return frozenset(
            tuple(row[k] for k in PROFILE_KEYS)
            for _, row in group.iterrows()
        )

    disease_profiles = (
        profile_df
        .groupby("ORPHAcode")
        .apply(build_profile_set)
        .rename("profile_set")
        .reset_index()
    )

    disease_names = df[["ORPHAcode", "DiseaseName"]].drop_duplicates("ORPHAcode")
    disease_profiles = disease_profiles.merge(disease_names, on="ORPHAcode", how="left")

    n_profiles = len(disease_profiles)

    # ------------------------------------------------------------------
    # 3. Build HPO sets H(d)
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 4. Classify targets vs donors
    # ------------------------------------------------------------------
    targets = diseases[diseases["n_hpo"] == 0].reset_index(drop=True)
    donors  = diseases[diseases["n_hpo"] >  0].reset_index(drop=True)

    n_targets = len(targets)
    n_donors  = len(donors)

    # ------------------------------------------------------------------
    # 5. Exact reaction-profile matching
    # ------------------------------------------------------------------
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
            results.append({
                "ORPHAcode_target":            d0_code,
                "DiseaseName_target":          d0_name,
                "n_exact_donors":              0,
                "ORPHAcode_donors":            "",
                "DiseaseName_donors":          "",
                "HPO_ID":                      "",
                "HPO_Label":                   "",
                "donor_count":                 0,
                "Supporting_ORPHAcode_donors": "",
                "candidate_type":              "no_match",
                "is_selected":                 False,
                "is_fully_supported":          False,
            })
            continue

        # ------------------------------------------------------------------
        # 6. HPO prioritization
        # ------------------------------------------------------------------
        hpo_donor_count: Dict[str, int]        = {}
        hpo_donor_codes: Dict[str, List[str]]  = {}

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

            is_selected        = ctype in ("dominant", "fully_supported")
            is_fully_supported = ctype == "fully_supported"

            results.append({
                "ORPHAcode_target":            d0_code,
                "DiseaseName_target":          d0_name,
                "n_exact_donors":              n_exact,
                "ORPHAcode_donors":            donor_codes_str,
                "DiseaseName_donors":          donor_names_str,
                "HPO_ID":                      hpo_id,
                "HPO_Label":                   hpo_label_map.get(hpo_id, ""),
                "donor_count":                 count,
                "Supporting_ORPHAcode_donors": "|".join(hpo_donor_codes[hpo_id]),
                "candidate_type":              ctype,
                "is_selected":                 is_selected,
                "is_fully_supported":          is_fully_supported,
            })

    # ------------------------------------------------------------------
    # 7. Assemble output
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(results)

    type_order = {
        "fully_supported": 0,
        "dominant":        1,
        "recurrent":       2,
        "single_donor":    3,
        "no_match":        4,
    }
    results_df["_type_order"] = results_df["candidate_type"].map(type_order)
    results_df = (
        results_df
        .sort_values(
            ["ORPHAcode_target", "_type_order", "donor_count"],
            ascending=[True, True, False],
        )
        .drop(columns=["_type_order"])
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------
    # 8. Compute summary statistics
    # ------------------------------------------------------------------
    n_matched  = results_df[results_df["candidate_type"] != "no_match"]["ORPHAcode_target"].nunique()
    n_no_match = results_df[results_df["candidate_type"] == "no_match"]["ORPHAcode_target"].nunique()

    # Per-type counts (target-level, not row-level)
    type_counts = {}
    for ct in ["fully_supported", "dominant", "recurrent", "single_donor", "no_match"]:
        candidates_ct = results_df[results_df["candidate_type"] == ct]
        type_counts[ct] = candidates_ct["ORPHAcode_target"].nunique() if len(candidates_ct) > 0 else 0

    # Candidate HPO terms per matched target
    matched_rows = results_df[results_df["candidate_type"] != "no_match"]
    candidates_per_target = matched_rows.groupby("ORPHAcode_target").size()
    mean_candidates = candidates_per_target.mean()
    median_candidates = candidates_per_target.median()
    min_candidates = candidates_per_target.min()
    max_candidates = candidates_per_target.max()

    # Total candidate HPO terms
    total_candidate_hpo = len(matched_rows)

    # Donors per matched target
    donors_per_target = matched_rows.groupby("ORPHAcode_target")["n_exact_donors"].first()
    mean_donors = donors_per_target.mean()

    stats = {
        "use_filter":                     use_cofactor_pathway_filter,
        "raw_rows":                       raw_rows,
        "raw_diseases":                   raw_diseases,
        "raw_diseases_with_hpo":          raw_diseases_with_hpo,
        "rows_after_filter":              rows_after_filter,
        "diseases_after_filter":          diseases_after_filter,
        "diseases_with_complete_profile": n_profiles,
        "n_targets":                      n_targets,
        "n_donors":                       n_donors,
        "n_matched":                      n_matched,
        "n_no_match":                     n_no_match,
        "type_counts":                    type_counts,
        "total_candidate_hpo":            total_candidate_hpo,
        "mean_candidates_per_target":     round(mean_candidates, 1),
        "median_candidates_per_target":   int(median_candidates),
        "min_candidates_per_target":      int(min_candidates),
        "max_candidates_per_target":      int(max_candidates),
        "mean_donors_per_matched_target": round(mean_donors, 1),
    }

    return results_df, stats


# ===================================================================
# MAIN: run both modes and print comparison
# ===================================================================
if __name__ == "__main__":
    print("=" * 72)
    print("Sensitivity analysis: cofactor/pathway filter ON vs OFF")
    print("=" * 72)

    all_results = {}
    all_stats   = {}

    for use_filter in [True, False]:
        label = "ON" if use_filter else "OFF"
        print(f"\n>>> Pipeline with filter = {label} ...\n")
        res, st = run_pipeline(use_filter)
        all_results[label] = res
        all_stats[label]   = st

    # ------------------------------------------------------------------
    # Print comparison table (plain text)
    # ------------------------------------------------------------------
    on  = all_stats["ON"]
    off = all_stats["OFF"]

    print("\n" + "=" * 72)
    print("COMPARISON SUMMARY")
    print("=" * 72)

    rows_comp = [
        ("Rows in input file",           f"{on['raw_rows']:,}",      f"{off['raw_rows']:,}",      "---"),
        ("Distinct diseases (raw)",      f"{on['raw_diseases']:,}",  f"{off['raw_diseases']:,}",  "---"),
        ("After Cof/Path filter",        f"{on['diseases_after_filter']:,}", f"{off['diseases_after_filter']:,}", "---"),
        ("With complete 6-field profile",f"{on['diseases_with_complete_profile']:,}", f"{off['diseases_with_complete_profile']:,}", f"{off['diseases_with_complete_profile'] - on['diseases_with_complete_profile']:+}"),
        ("",                             "",                          "",                          ""),
        ("Targets (H=empty)",            f"{on['n_targets']:,}",     f"{off['n_targets']:,}",     f"{off['n_targets'] - on['n_targets']:+}"),
        ("Donors (|H|>0)",               f"{on['n_donors']:,}",      f"{off['n_donors']:,}",      f"{off['n_donors'] - on['n_donors']:+}"),
        ("",                             "",                          "",                          ""),
        ("Matched targets (exact match)",f"{on['n_matched']:,}",     f"{off['n_matched']:,}",     f"{off['n_matched'] - on['n_matched']:+}"),
        ("No-match targets",             f"{on['n_no_match']:,}",    f"{off['n_no_match']:,}",    f"{off['n_no_match'] - on['n_no_match']:+}"),
        ("",                             "",                          "",                          ""),
    ]

    # Add per-type counts
    for ct in ["fully_supported", "dominant", "recurrent", "single_donor", "no_match"]:
        rows_comp.append(
            (f"  Targets with {ct} candidates",
             f"{on['type_counts'][ct]:,}",
             f"{off['type_counts'][ct]:,}",
             f"{off['type_counts'][ct] - on['type_counts'][ct]:+}")
        )

    rows_comp += [
        ("",                             "",                          "",                          ""),
        ("Total candidate HPO terms",    f"{on['total_candidate_hpo']:,}", f"{off['total_candidate_hpo']:,}", f"{off['total_candidate_hpo'] - on['total_candidate_hpo']:+}"),
        ("Mean candidates / matched target", f"{on['mean_candidates_per_target']}", f"{off['mean_candidates_per_target']}", "---"),
        ("Median candidates / matched target", f"{on['median_candidates_per_target']}", f"{off['median_candidates_per_target']}", "---"),
        ("Min / Max candidates",         f"{on['min_candidates_per_target']}/{on['max_candidates_per_target']}", f"{off['min_candidates_per_target']}/{off['max_candidates_per_target']}", "---"),
        ("Mean donors / matched target", f"{on['mean_donors_per_matched_target']}", f"{off['mean_donors_per_matched_target']}", "---"),
    ]

    header = f"{'Metric':<55} {'Filter ON':>12} {'Filter OFF':>12} {'Delta':>8}"
    sep = "-" * len(header)
    print(f"\n{header}")
    print(sep)
    for label, v_on, v_off, delta in rows_comp:
        if not label:
            print()
        else:
            print(f"{label:<55} {v_on:>12} {v_off:>12} {delta:>8}")

    # ------------------------------------------------------------------
    # LaTeX table for the paper
    # ------------------------------------------------------------------
    print("\n" + "=" * 72)
    print("LATEX TABLE (copy-paste into paper)")
    print("=" * 72)

    # Helper to format delta with sign
    def fmt_delta(val_on, val_off):
        d = val_off - val_on
        if d > 0:
            return f"$+{d}$"
        elif d < 0:
            return f"${d}$"
        else:
            return "$0$"

    print("""
\\begin{table}[t]
\\centering
\\caption{Sensitivity analysis of the cofactor/pathway filter.
         Filter~ON requires each disease to have at least one annotated
         cofactor and at least one pathway; Filter~OFF omits this
         requirement.  The 6-field reaction profile
         (Gene, Entry, EC, Rhea\textunderscore ID, Cofactor, Pathway)
         and complete-profile requirement ($\\Phi_j(d)$ non-null for all
         six fields) are applied in both settings.}
\\label{tab:sensitivity}
\\small
\\begin{tabular}{lrrr}
\\toprule
\\textbf{Metric} & \\textbf{Filter ON} & \\textbf{Filter OFF} & $\\bm{\\Delta}$ \\\\
\\midrule""")

    latex_rows = [
        ("Diseases with complete 6-field profile",
         f"{on['diseases_with_complete_profile']:,}",
         f"{off['diseases_with_complete_profile']:,}",
         fmt_delta(on['diseases_with_complete_profile'], off['diseases_with_complete_profile'])),
        ("Target diseases ($H(d)=\\\\emptyset$)",
         f"{on['n_targets']:,}",
         f"{off['n_targets']:,}",
         fmt_delta(on['n_targets'], off['n_targets'])),
        ("Donor diseases ($|H(d)|>0$)",
         f"{on['n_donors']:,}",
         f"{off['n_donors']:,}",
         fmt_delta(on['n_donors'], off['n_donors'])),
        ("\\midrule"),
        ("Matched targets (exact reaction-profile match)",
         f"{on['n_matched']:,}",
         f"{off['n_matched']:,}",
         fmt_delta(on['n_matched'], off['n_matched'])),
        ("No-match targets",
         f"{on['n_no_match']:,}",
         f"{off['n_no_match']:,}",
         fmt_delta(on['n_no_match'], off['n_no_match'])),
        ("\\midrule"),
        ("Total candidate HPO terms",
         f"{on['total_candidate_hpo']:,}",
         f"{off['total_candidate_hpo']:,}",
         fmt_delta(on['total_candidate_hpo'], off['total_candidate_hpo'])),
    ]

    for ct_label, ct_key in [("fully\\_supported", "fully_supported"),
                              ("dominant", "dominant"),
                              ("recurrent", "recurrent"),
                              ("single\\_donor", "single_donor")]:
        latex_rows.append(
            (f"~~Targets with {ct_label} candidates",
             f"{on['type_counts'][ct_key]:,}",
             f"{off['type_counts'][ct_key]:,}",
             fmt_delta(on['type_counts'][ct_key], off['type_counts'][ct_key]))
        )

    latex_rows += [
        ("\\midrule"),
        ("Mean candidates per matched target",
         str(on['mean_candidates_per_target']),
         str(off['mean_candidates_per_target']),
         "---"),
        ("Median candidates per matched target",
         str(on['median_candidates_per_target']),
         str(off['median_candidates_per_target']),
         "---"),
        ("Mean donor diseases per matched target",
         str(on['mean_donors_per_matched_target']),
         str(off['mean_donors_per_matched_target']),
         "---"),
    ]

    for parts in latex_rows:
        if isinstance(parts, str) and parts.startswith("\\midrule"):
            print("\\midrule")
        elif isinstance(parts, tuple) and parts[0].startswith("\\midrule"):
            print("\\midrule")
        else:
            lbl, v_on, v_off, d = parts
            print(f"{lbl} & {v_on} & {v_off} & {d} \\\\")

    print("""\\bottomrule
\\end{tabular}
\\end{table}""")

    # ------------------------------------------------------------------
    # Also save both outputs to CSV for inspection
    # ------------------------------------------------------------------
    out_on  = r"C:\Users\biofa\Desktop\DATAEXTRACTION\hpo_prioritization_filter_ON.csv"
    out_off = r"C:\Users\biofa\Desktop\DATAEXTRACTION\hpo_prioritization_filter_OFF.csv"
    all_results["ON"].to_csv(out_on,  sep=";", index=False, encoding="utf-8-sig")
    all_results["OFF"].to_csv(out_off, sep=";", index=False, encoding="utf-8-sig")
    print(f"\nDetailed results saved to:")
    print(f"  Filter ON : {out_on}")
    print(f"  Filter OFF: {out_off}")
    print("\nDone.")
