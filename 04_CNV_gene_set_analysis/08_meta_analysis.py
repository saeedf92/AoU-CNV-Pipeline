"""
Cross-Ancestry Meta-Analysis Pipeline: Rare CNV Gene Set Burden
================================================================
Performs weighted Stouffer's Z meta-analysis across EUR, AFR, and AMR
ancestries and computes Cochran's Q / I² heterogeneity statistics.

HOW TO CONFIGURE YOUR OWN ANALYSIS
------------------------------------
(IMPORTANT) Edit the USER CONFIGURATION section below. At minimum you need to:

  1. Set N_SAMPLES for your cohort sizes.
  2. Define one or more AnalysisConfig entries in ANALYSES.
  3. Supply your target gene sets either as:
       a) A plain-text file (one gene set name per line) → set targets_file
       b) A Python list passed directly                 → set targets
       c) Leave both None to include ALL gene sets found in the result files.

Usage:
  python meta_analysis_pipeline.py
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from scipy.stats import norm, chi2
from statsmodels.stats.multitest import multipletests


# =============================================================================
# USER CONFIGURATION  ← edit this section
# =============================================================================

# Sample sizes per ancestry label (must match your folder names)
N_SAMPLES: dict[str, int] = {
    "EUR": 34868,
    "AFR": 12720,
    "AMR": 12353,
}

# Ancestry labels (must match subfolder names: EUR/, AFR/, AMR/)
ANCESTRIES: list[str] = ["EUR", "AFR", "AMR"]


@dataclass
class AnalysisConfig:
    """
    Describes a single meta-analysis run.

    Parameters
    ----------
    label : str
        Short human-readable name shown in output headers (e.g. "NDD PGC53").

    model_path : str
        Path *inside* each ancestry folder to the directory holding result files
        (e.g. "ndd_pgc_rare_cnv_count/pgc_covs_age_sex/firth").

    output_path : str
        Where to write the combined TSV result.

    input_file : str or None
        If your results are pre-merged into a single file per ancestry (e.g.
        "PGC53_6METRICS.tsv"), set this to that filename.
        Leave None to load individual per-gene-set TSV files from model_path.

    targets : list[str] or None
        Gene set names to include. Matched fuzzily against filenames or the
        Match_ID column. Pass an empty list or None to include everything found.

    targets_file : str or None
        Path to a plain-text file with one gene set name per line.
        Takes precedence over `targets` if both are supplied.

    use_10kb : bool
        Convenience toggle — if True, appends "_10kb" to the size tag in
        the summary header. Does not change any paths; those must be set
        explicitly in model_path / output_path.

    report_metric_overlap : bool
        Whether to print the per-CNV-metric FDR overlap table (useful for
        analyses that include dup/del count and burden metrics).
    """
    label: str
    model_path: str
    output_path: str
    input_file: str | None = None
    targets: list[str] = field(default_factory=list)
    targets_file: str | None = None
    use_10kb: bool = False
    report_metric_overlap: bool = False


# ---------------------------------------------------------------------------
# Define your analyses here.  Add, remove, or modify entries as needed.
# ---------------------------------------------------------------------------
ANALYSES: list[AnalysisConfig] = [

    # -- NDD: 53 PGC core gene sets (pre-merged summary file per ancestry) --
    AnalysisConfig(
        label="NDD PGC53",
        model_path="ndd_pgc_rare_cnv_count/pgc_covs_age_sex/firth",
        output_path="meta_analysis/ndd/META_STOUFFERS_PGC53_NDD_ALL.tsv",
        input_file="PGC53_6METRICS.tsv",
        # Supply targets as a list. # For example: 
        targets=[  
            "PhHs_NervSys_All", "Neurof_GoNeuronProj", "Neurof_UnionStringent", "Netw_full",
            "PhHs_NervSys_ADX", "scRNA_Expressed_PgS", "scRNA_Expressed_ExDp2", "Kirov_ARC",
            "BspanVH_lg2rpkm4.74", "Neurof_GoSynaptic", "PSD_BayesGrant_fullset", "BspanHM_lg2rpkm3.21",
            "PhMm_NervSystem_all", "scRNA_Expressed_vRG", "gnomAD_oe_lof_upper_0.35", "scRNA_Expressed_ExM",
            "scRNA_Expressed_PgG2M", "scRNA_Expressed_ExM_U", "scRNA_Expressed_oRG", "Kirov_NMDAR",
            "scRNA_Expressed_OPC", "BspanVHM_PreNat", "BspanVHM_PstNat", "scRNA_Expressed_ExDp1",
            "scRNA_Expressed_Per", "Neurof_GoNervTransm", "Neurof_GoNervSysDev", "Neurof_GoNeuronBody",
            "BspanLA_lg2rpkm.MIN", "scRNA_Expressed_End", "scRNA_Expressed_Mic", "scRNA_Expressed_IP",
            "PhMm_NeuroUnion_all", "Neurof_UnionInclusive", "scRNA_Expressed_InCGE", "scRNA_Expressed_ExN",
            "PhMm_Aggr_EndoExocrRepr_all", "scRNA_Expressed_InMGE", "PhMm_Aggr_Sensory_all",
            "PhMm_Aggr_CardvascMuscle_all", "FMR1_Targets_Darnell", "FMR1_Targets_Ascano",
            "PhMm_NeuroBehav_all", "PhMm_Aggr_IntegAdipPigm_all", "Neurof_GoNervSysDev_CNS",
            "BspanML_lg2rpkm0.93", "Neurof_PathwaysAxonG", "PhHs_MindFun_All",
            "PhMm_Aggr_HematoImmune_all", "Neurof_KeggSynaptic", "PhHs_MindFun_ADX",
            "PhMm_Aggr_DigestHepato_all", "PhMm_Aggr_SkeCranioLimbs_all",
        ],
        # Or point to a text file instead (one name per line):
        # targets_file="my_gene_sets.txt",
    ),

    # -- NDD 10kb (same targets, different model path) --
    AnalysisConfig(
        label="NDD PGC53 10KB",
        model_path="ndd_pgc_rare_10kb_cnv_count/pgc_covs_age_sex/firth",
        output_path="meta_analysis/ndd/META_STOUFFERS_PGC53_NDD_10KB.tsv",
        input_file="PGC53_6METRICS.tsv",
        use_10kb=True,
        # Leaving targets empty → includes all gene sets found in the file
    ),

    # -- PTSD GWAS gene sets (individual TSV files per gene set) --
    AnalysisConfig(
        label="PTSD GWAS",
        model_path="pgc_ptsd_gwas_rare_cnv_count/pgc_covs_age_sex/firth",
        output_path="meta_analysis/pgc_ptsd/META_STOUFFERS_PGC_PTSD_ALL.tsv",
        # No input_file → loads individual TSVs from the folder
        # No targets    → includes all files found
    ),

    # -- Abnormal phenotype gene sets (pre-merged summary file) --
    AnalysisConfig(
        label="Abnormal 146",
        model_path="abnormal_gene_sets_rare_cnv_count/pgc_covs_age_sex/firth",
        output_path="meta_analysis/abnormal/META_STOUFFERS_ABNORMAL146_ALL.tsv",
        input_file="ABNORMAL146_6METRICS.tsv",
        report_metric_overlap=True,
    ),

]


# =============================================================================
# INTERNAL UTILITIES  (no editing needed below this line)
# =============================================================================

def _load_targets(cfg: AnalysisConfig) -> list[str]:
    """
    Resolves the target gene set list from an AnalysisConfig.
    Priority: targets_file > targets list > empty (all sets).
    """
    if cfg.targets_file:
        path = cfg.targets_file
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"targets_file not found: {path!r}\n"
                "Create a plain-text file with one gene set name per line."
            )
        with open(path) as fh:
            names = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        print(f"  Loaded {len(names)} target gene sets from {path}")
        return names
    return list(cfg.targets)  # may be empty → no filtering


def _clean_key(text: str) -> str:
    """Normalises a string for fuzzy filename / ID matching."""
    if not isinstance(text, str):
        return ""
    for token in ["GW_EUR-", "GW_AFR-", "GW_AMR-", "_burden_summary_results",
                  "_summary_results", "_10kb_summary_results", ".tsv"]:
        text = text.replace(token, "")
    return "".join(c for c in text.lower() if c.isalnum())


def _clean_key_series(series: pd.Series) -> pd.Series:
    """Vectorised version of _clean_key for DataFrame columns."""
    return series.astype(str).str.strip().str.lower().str.replace(r"[^a-z0-9]+", "", regex=True)


def _load_tsv(path: str, match_id: str | None = None) -> pd.DataFrame | None:
    """Reads a result TSV, normalises column names, returns tidy subset or None."""
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        print(f"  [WARN] Could not read {path}: {exc}")
        return None

    # Accept alternative column names for the gene-set identifier
    for alt in ("PGC_NAME", "Gene_Set_Label", "gene_set", "gene_set_name"):
        if alt in df.columns and "Match_ID" not in df.columns:
            df = df.rename(columns={alt: "Match_ID"})

    if "Match_ID" not in df.columns:
        if match_id is not None:
            df["Match_ID"] = match_id
        else:
            print(f"  [WARN] {path}: no Match_ID column and no match_id supplied — skipping.")
            return None

    required = ["Match_ID", "Gene_Set_Metric", "Beta", "Std_Error", "P_Value"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"  [WARN] {path} missing columns: {missing} — skipping.")
        return None

    return df[required].copy()


def _clean_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Coerces Beta/SE/P to numeric, drops non-finite rows, requires SE > 0."""
    for col in ("Beta", "Std_Error", "P_Value"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Beta", "Std_Error", "P_Value"])
    df = df[np.isfinite(df[["Beta", "Std_Error", "P_Value"]]).all(axis=1)]
    return df[df["Std_Error"] > 0].copy()


# ---------------------------------------------------------------------------
# Data loading strategies
# ---------------------------------------------------------------------------

def _load_from_single_file(
    cfg: AnalysisConfig,
    targets: list[str],
) -> dict[str, pd.DataFrame]:
    """Loads one pre-merged TSV per ancestry and optionally filters to targets."""
    target_set = set(targets)
    pop_dfs: dict[str, pd.DataFrame] = {}

    for anc in ANCESTRIES:
        path = os.path.join(anc, cfg.model_path, cfg.input_file)
        print(f"  {anc}: {path}  exists={os.path.exists(path)}")
        df = _load_tsv(path)
        if df is None:
            continue
        df = _clean_numeric(df)
        if target_set:
            df = df[df["Match_ID"].isin(target_set)]
            if df.empty:
                print(f"  [WARN] {anc}: no rows matched the target list — check gene set names.")
                continue
        # Fuzzy join keys (needed when the Abnormal module uses Match_Key/Metric_Key)
        df["Match_Key"] = _clean_key_series(df["Match_ID"])
        df["Metric_Key"] = _clean_key_series(df["Gene_Set_Metric"])
        pop_dfs[anc] = df

    return pop_dfs


def _load_from_individual_files(
    cfg: AnalysisConfig,
    targets: list[str],
) -> dict[str, pd.DataFrame]:
    """Loads individual per-gene-set TSV files from each ancestry folder."""
    target_lookup = {_clean_key(t): t for t in targets} if targets else {}
    pop_dfs: dict[str, pd.DataFrame] = {}

    for anc in ANCESTRIES:
        folder = os.path.join(anc, cfg.model_path)
        all_files = [
            f for f in glob.glob(os.path.join(folder, "*.tsv"))
            if "summary_results" in f.lower() and "FINAL" not in os.path.basename(f).upper()
        ]
        rows = []
        for f_path in all_files:
            # Derive a clean gene-set name from the filename
            name = (
                os.path.basename(f_path)
                .replace(f"GW_{anc}_", "")
                .replace("_summary_results.tsv", "")
                .replace("_10kb_summary_results.tsv", "")
            )
            # If a target list was given, skip files that don't match
            if target_lookup:
                canonical = target_lookup.get(_clean_key(name))
                if canonical is None:
                    continue
                name = canonical

            df = _load_tsv(f_path, match_id=name)
            if df is not None:
                rows.append(df)

        if rows:
            combined = pd.concat(rows, ignore_index=True)
            combined = _clean_numeric(combined)
            combined["Match_Key"] = _clean_key_series(combined["Match_ID"])
            combined["Metric_Key"] = _clean_key_series(combined["Gene_Set_Metric"])
            pop_dfs[anc] = combined

    return pop_dfs


# ---------------------------------------------------------------------------
# Merging and statistics
# ---------------------------------------------------------------------------

def _merge_ancestries(pop_dfs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Inner-joins ancestry DataFrames on fuzzy Match_Key + Metric_Key."""
    # Rename EUR columns, preserving the readable ID for output
    eur = pop_dfs["EUR"].rename(
        columns={
            "Beta": "B_EUR", "Std_Error": "SE_EUR", "P_Value": "P_EUR",
            "Match_ID": "Match_ID_EUR", "Gene_Set_Metric": "Gene_Set_Metric_EUR",
        }
    )
    combined = eur

    for anc in [a for a in ANCESTRIES if a != "EUR"]:
        if anc not in pop_dfs:
            continue
        other = pop_dfs[anc].rename(
            columns={
                "Beta": f"B_{anc}", "Std_Error": f"SE_{anc}", "P_Value": f"P_{anc}",
                "Match_ID": f"Match_ID_{anc}", "Gene_Set_Metric": f"Gene_Set_Metric_{anc}",
            }
        )
        before = len(combined)
        combined = pd.merge(combined, other, on=["Match_Key", "Metric_Key"], how="inner")
        print(f"  Merged {anc}: {before} → {len(combined)} tests")

    # Restore human-readable labels from EUR
    combined["Match_ID"] = combined["Match_ID_EUR"]
    combined["Gene_Set_Metric"] = combined["Gene_Set_Metric_EUR"]
    return combined


def _stouffers_meta(combined: pd.DataFrame, active_pops: list[str]) -> pd.DataFrame:
    """
    Computes Stouffer's Z, IVW beta, Cochran's Q, I², P_het, FDR, and direction.
    Returns a copy of the input with new columns appended.
    """
    df = combined.copy()
    total_n = sum(N_SAMPLES[p] for p in active_pops)

    # Weighted Stouffer's Z (sample-size weights)
    weighted_z = sum(
        (df[f"B_{p}"] / df[f"SE_{p}"]) * np.sqrt(N_SAMPLES[p])
        for p in active_pops
    )
    df["Meta_Z"] = weighted_z / np.sqrt(total_n)
    df["Meta_P_Weighted"] = 2 * (1 - norm.cdf(df["Meta_Z"].abs()))

    # Inverse-Variance Weighted (IVW) beta
    sum_w  = sum(1 / df[f"SE_{p}"] ** 2 for p in active_pops)
    sum_wb = sum(df[f"B_{p}"] / df[f"SE_{p}"] ** 2 for p in active_pops)
    b_ivw  = sum_wb / sum_w

    df["Meta_Beta_IVW"] = b_ivw
    df["Meta_SE_IVW"]   = np.sqrt(1 / sum_w)
    df["Meta_Z_IVW"]    = df["Meta_Beta_IVW"] / df["Meta_SE_IVW"]
    df["Meta_P_IVW"]    = 2 * (1 - norm.cdf(df["Meta_Z_IVW"].abs()))

    # Cochran's Q and I²
    k = len(active_pops)
    q = sum((1 / df[f"SE_{p}"] ** 2) * (df[f"B_{p}"] - b_ivw) ** 2 for p in active_pops)
    df["Cochrans_Q"] = q
    df["P_het"]      = 1 - chi2.cdf(q, df=k - 1)
    df["I_squared"]  = np.where(q > 0, np.maximum(0, (q - (k - 1)) / q) * 100, 0)

    # FDR correction (Benjamini-Hochberg)
    df["Meta_Q_Weighted"] = multipletests(df["Meta_P_Weighted"], method="fdr_bh")[1]
    df["Meta_Q_IVW"]      = multipletests(df["Meta_P_IVW"],      method="fdr_bh")[1]

    # Direction consistency string  e.g. "++−"
    def _direction(row: pd.Series) -> str:
        return "".join(
            "+" if row[f"B_{p}"] > 0 else ("-" if row[f"B_{p}"] < 0 else "0")
            for p in active_pops
        )
    df["Direction"] = df.apply(_direction, axis=1)

    return df


def _save(df: pd.DataFrame, output_path: str) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    df.sort_values("Meta_P_Weighted").to_csv(output_path, sep="\t", index=False)
    print(f"  Results saved → {output_path}")


# =============================================================================
# REPORTING
# =============================================================================

def _print_summary(df: pd.DataFrame, label: str) -> None:
    sig_p = df[df["Meta_P_Weighted"] < 0.05]
    sig_q = df[df["Meta_Q_Weighted"] < 0.05]

    sep = "=" * 80
    print(f"\n{sep}")
    print(f"  META-ANALYSIS SUMMARY — {label}")
    print(sep)
    print(f"  Unique gene sets:         {df['Match_ID'].nunique()}")
    print(f"  Total tests:              {len(df)}")
    print(f"  {'─' * 76}")
    print(f"  Nominal  (P < 0.05):      {len(sig_p)} tests  |  {sig_p['Match_ID'].nunique()} unique sets")
    print(f"  FDR      (Q < 0.05):      {len(sig_q)} tests  |  {sig_q['Match_ID'].nunique()} unique sets")
    print(f"  High heterogeneity I²>50: {(df['I_squared'] > 50).sum()} tests")
    print(sep)

    if not sig_p.empty:
        cols = ["Match_ID", "Gene_Set_Metric", "Meta_P_Weighted",
                "Meta_Q_Weighted", "I_squared", "P_het"]
        available = [c for c in cols if c in df.columns]
        print("\n  TOP 10 NOMINAL HITS")
        print("  " + "─" * 76)
        print(sig_p[available].sort_values("Meta_P_Weighted").head(10).to_string(index=False))
    else:
        print("\n  No nominally significant hits found.")


def _print_direction_summary(df: pd.DataFrame, active_pops: list[str]) -> None:
    all_pos = all((df[f"B_{p}"] > 0).all() for p in active_pops)  # noqa — used below
    consistent_mask = (
        df[[f"B_{p}" for p in active_pops]].gt(0).all(axis=1) |
        df[[f"B_{p}" for p in active_pops]].lt(0).all(axis=1)
    )
    n_consistent = consistent_mask.sum()
    pct = 100 * n_consistent / len(df) if len(df) else 0

    robust = df[consistent_mask & (df["Meta_P_Weighted"] < 0.05)]
    print(f"\n  Directional consistency:  {n_consistent}/{len(df)} ({pct:.1f}%) tests")
    print(f"  Robust hits (P<0.05 + consistent direction): {len(robust)}")

    if not robust.empty:
        cols = (["Match_ID", "Gene_Set_Metric"]
                + [f"B_{p}" for p in active_pops]
                + ["Meta_P_Weighted", "Meta_Q_Weighted", "I_squared", "P_het"])
        available = [c for c in cols if c in df.columns]
        print("\n  TOP 5 ROBUST HITS")
        print("  " + "─" * 76)
        print(
            robust[available].sort_values("Meta_P_Weighted").head(5)
            .to_string(index=False, float_format=lambda x: f"{x:.4f}")
        )


def _print_metric_overlap(df: pd.DataFrame) -> None:
    """Per-CNV-metric FDR overlap table (count vs. burden, dup vs. del)."""
    metrics = {
        "dup_count":  "genome_wide_dup_count",
        "dup_burden": "genome_wide_dup_burden_mb",
        "del_count":  "genome_wide_del_count",
        "del_burden": "genome_wide_del_burden_mb",
    }
    print(f"\n  {'─' * 76}")
    print("  CNV METRIC FDR OVERLAP (Q < 0.05)")
    print(f"  {'─' * 76}")

    for framework, q_col in [("STOUFFER", "Meta_Q_Weighted"), ("IVW", "Meta_Q_IVW")]:
        sets = {
            k: set(df.loc[(df[q_col] < 0.05) & (df["Gene_Set_Metric"] == v), "Match_ID"])
            for k, v in metrics.items()
        }
        print(f"\n  [{framework}]")
        print(f"    DUP  count={len(sets['dup_count']):<4}  burden={len(sets['dup_burden']):<4}  "
              f"overlap={len(sets['dup_count'] & sets['dup_burden'])}")
        print(f"    DEL  count={len(sets['del_count']):<4}  burden={len(sets['del_burden']):<4}  "
              f"overlap={len(sets['del_count'] & sets['del_burden'])}")


# =============================================================================
# MAIN RUNNER
# =============================================================================

def run_analysis(cfg: AnalysisConfig) -> None:
    """Executes a single meta-analysis run described by an AnalysisConfig."""
    print(f"\n{'─' * 80}")
    print(f"  Running: {cfg.label}")
    print(f"{'─' * 80}")

    # Resolve target gene sets
    targets = _load_targets(cfg)
    if targets:
        print(f"  Target gene sets: {len(targets)}")
    else:
        print("  Target gene sets: ALL (no filter applied)")

    # Load data per ancestry
    if cfg.input_file:
        pop_dfs = _load_from_single_file(cfg, targets)
    else:
        pop_dfs = _load_from_individual_files(cfg, targets)

    if "EUR" not in pop_dfs:
        print(f"  [ERROR] EUR data not found for {cfg.label!r}. Skipping.")
        return

    # Merge and compute
    combined = _merge_ancestries(pop_dfs)
    if combined.empty:
        print(f"  [ERROR] No overlapping tests after merge for {cfg.label!r}. Skipping.")
        return

    active_pops = [a for a in ANCESTRIES if f"B_{a}" in combined.columns]
    print(f"  Meta-analyzing {len(combined)} overlapping tests across {active_pops}")

    result = _stouffers_meta(combined, active_pops)
    _save(result, cfg.output_path)
    _print_summary(result, cfg.label)
    _print_direction_summary(result, active_pops)

    if cfg.report_metric_overlap:
        _print_metric_overlap(result)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-ancestry CNV burden meta-analysis pipeline."
    )
    parser.add_argument(
        "--analyses", nargs="*", metavar="LABEL",
        help="Run only the named analysis labels (space-separated). "
             "Omit to run all analyses defined in ANALYSES.",
    )
    args = parser.parse_args()

    selected = ANALYSES
    if args.analyses:
        label_set = set(args.analyses)
        selected = [a for a in ANALYSES if a.label in label_set]
        if not selected:
            print(f"[ERROR] No analyses matched labels: {args.analyses}")
            available = [a.label for a in ANALYSES]
            print(f"Available: {available}")
            return

    print("\n" + "=" * 80)
    print("  CROSS-ANCESTRY CNV BURDEN META-ANALYSIS")
    print("=" * 80)

    for cfg in selected:
        run_analysis(cfg)

    print(f"\n{'=' * 80}")
    print("  Done.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()