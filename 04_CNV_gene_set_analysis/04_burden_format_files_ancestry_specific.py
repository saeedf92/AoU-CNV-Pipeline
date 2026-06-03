"""
CNV Burden Analysis Pipeline
=============================
Universal pipeline for formatting and filtering CNV annotations per ancestry,
and extracting rare CNV IDs overlapping gene sets (AUD, NDD, PGC GWAS, Abnormal,
Housekeeping/negative control).

Usage
-----
Configure the settings block at the bottom of this file, then run:
    python cnv_burden_pipeline.py

Dependencies: pandas, matplotlib, os, re
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# CONFIGURATION
# =============================================================================

WORKSPACE_BUCKET = os.getenv("WORKSPACE_BUCKET", "")

ANNOTATION_FILES = {
    "EUR": "./cnv_vcf_plink/AoU_srWGS_SV.v8_CNVs_EUR_final_master_annotation.tsv",
    "AFR": "./cnv_vcf_plink/AoU_srWGS_SV.v8_CNVs_AFR_final_master_annotation.tsv",
    "AMR": "./cnv_vcf_plink/AoU_srWGS_SV.v8_CNVs_AMR_final_master_annotation.tsv",
}

GENE_SET_FILES = {
    "aud":          "GENE-SETS/aud_gene_sets.txt",
    "ndd":          "GENE-SETS/pgc_ptsd_cnv_ndd_gene_sets.txt",
    "housekeeping": "GENE-SETS/gene_housekeeping_set_neg_control.txt",
    "pgc_gwas": [
        "GENE-SETS/PGC_GWAS_GENES/brain_tissue_eQTL_smr_suppl_tbl_17.txt",
        "GENE-SETS/PGC_GWAS_GENES/magma_suppl_tbl_13.txt",
        "GENE-SETS/PGC_GWAS_GENES/brain_tissue_sc_dlPFC_suppl_tbl_19.txt",
        "GENE-SETS/PGC_GWAS_GENES/protein_coding_suppl_tbl_9.txt",
        "GENE-SETS/PGC_GWAS_GENES/brain_tissue_twas_suppl_tbl_16.txt",
    ],
    "abnormal_dir": "GENE-SETS/ListGeneSetsN",
}

# Columns in annotation files
CNV_ID_COL   = "id"
GENE_COL     = "collapsed_genes"
POP_AF_COL   = "pop_af"
LENGTH_COL   = "total_cnv_length_bp"
OVERLAP_COL  = "cnv_gene_overlap_length"
CHR_COL      = "chr"

AUTOSOMES    = [str(i) for i in range(1, 23)]
RARE_AF_CUTOFF = 0.01
SIZE_10KB      = 10_000

OUTPUT_ROOT  = "GENE-SETS"
PLOT_DIR     = os.path.join(OUTPUT_ROOT, "plots")


# =============================================================================
# HELPERS
# =============================================================================

def load_annotations(annotation_files: dict) -> dict:
    """Load per-ancestry annotation TSVs into a dict of DataFrames."""
    return {
        anc: pd.read_csv(path, sep="\t")
        for anc, path in annotation_files.items()
    }


def filter_rare_del_dup(df: pd.DataFrame, min_size_bp: int = 0) -> pd.DataFrame:
    """
    Return autosomal rare (pop_af < RARE_AF_CUTOFF) DEL/DUP CNVs.
    Optionally restrict to CNVs >= min_size_bp.
    """
    df = df.copy()
    df[POP_AF_COL]  = pd.to_numeric(df[POP_AF_COL],  errors="coerce")
    df[LENGTH_COL]  = pd.to_numeric(df[LENGTH_COL],  errors="coerce")
    df["_chr"]      = df[CHR_COL].astype(str).str.replace("chr", "", case=False)

    mask = (
        df["_chr"].isin(AUTOSOMES)
        & (df[POP_AF_COL] < RARE_AF_CUTOFF)
        & df[CNV_ID_COL].astype(str).str.contains("DEL|DUP", na=False)
    )
    if min_size_bp > 0:
        mask &= df[LENGTH_COL] >= min_size_bp

    return df[mask].drop(columns=["_chr"])


def expand_genes(df: pd.DataFrame, gene_col: str = GENE_COL) -> pd.DataFrame:
    """Explode semicolon-delimited gene names into one row per gene."""
    df = df.copy()
    df[gene_col] = df[gene_col].fillna(".").astype(str)
    df = df[df[gene_col] != "."]
    df[gene_col] = df[gene_col].str.split(";")
    df = df.explode(gene_col)
    df[gene_col] = df[gene_col].str.strip().str.upper()
    return df[df[gene_col] != ""]


def safe_name(x: str) -> str:
    """Convert an arbitrary string to a filesystem-safe identifier."""
    x = str(x).replace(".txt", "")
    x = re.sub(r"\s+", "_", x)
    x = re.sub(r"[^\w.-]", "_", x)
    return x


def save_cnv_ids(series: pd.Series, path: str) -> None:
    """Deduplicate, sort, and write a CNV ID list (no header)."""
    series.dropna().drop_duplicates().sort_values().to_csv(
        path, index=False, header=False
    )


def get_overlapping(cnv_gene_df: pd.DataFrame, gene_set: set) -> tuple:
    """
    Return (overlapping_cnv_ids, overlapping_gene_names) for a given gene set.
    """
    overlap = cnv_gene_df[cnv_gene_df[GENE_COL].isin(gene_set)]
    cnv_ids = overlap[CNV_ID_COL].dropna().drop_duplicates()
    genes   = overlap[GENE_COL].dropna().drop_duplicates()
    return cnv_ids, genes


# =============================================================================
# MODULE 1: RARE ANNOTATION EXPORT
# =============================================================================

def export_rare_annotations(annotations: dict, output_dir: str) -> None:
    """Save per-ancestry rare DEL/DUP annotation TSVs."""
    os.makedirs(output_dir, exist_ok=True)
    for anc, df in annotations.items():
        out = filter_rare_del_dup(df)
        path = os.path.join(output_dir, f"AoU_srWGS_SV.v8_{anc}_RARE_ANNOTATIONS.txt")
        out.to_csv(path, sep="\t", index=False)
        print(f"[{anc}] Rare DEL/DUP rows: {len(out):,}  →  {path}")


# =============================================================================
# MODULE 2: SIZE STATISTICS
# =============================================================================

def print_size_stats(annotations: dict) -> pd.DataFrame:
    """Print and return a summary table of CNV length distributions."""
    rows = []
    for anc, df in annotations.items():
        rare = filter_rare_del_dup(df)
        lt10kb  = rare[rare[LENGTH_COL] < 10_000]
        lt100bp = rare[rare[LENGTH_COL] < 100]
        rows.append({
            "Ancestry":              anc,
            "N_rare_DEL_DUP":        len(rare),
            "N_lt_10kb":             len(lt10kb),
            "Pct_lt_10kb":           round(len(lt10kb) / len(rare) * 100, 2) if len(rare) else 0,
            "N_lt_100bp":            len(lt100bp),
            "Pct_lt_100bp":          round(len(lt100bp) / len(rare) * 100, 2) if len(rare) else 0,
            "Median_length_bp":      rare[LENGTH_COL].median(),
            "Mean_length_bp":        round(rare[LENGTH_COL].mean(), 2),
            "Max_length_bp":         rare[LENGTH_COL].max(),
        })
    stats_df = pd.DataFrame(rows)
    print(stats_df.to_string(index=False))
    return stats_df


# =============================================================================
# MODULE 3: EXTRACT RARE CNV IDs BY SIZE CATEGORY
# =============================================================================

def export_cnv_id_lists(annotations: dict, output_dir: str) -> None:
    """
    Export per-ancestry CNV ID lists for:
      - all rare DEL/DUPs
      - rare >=10kb DEL/DUPs
      - rare genic DEL/DUPs
    """
    os.makedirs(output_dir, exist_ok=True)
    for anc, df in annotations.items():
        rare      = filter_rare_del_dup(df)
        rare_10kb = filter_rare_del_dup(df, min_size_bp=SIZE_10KB)

        # Genic filter
        rare[OVERLAP_COL] = pd.to_numeric(rare.get(OVERLAP_COL, pd.Series(dtype=float)), errors="coerce")
        genic_mask = (
            (rare[GENE_COL].notna() & (rare[GENE_COL] != "."))
            | (rare[OVERLAP_COL] >= 1)
        )
        genic = rare[genic_mask]

        save_cnv_ids(rare[CNV_ID_COL],      os.path.join(output_dir, f"{anc}_rare_cnvs.txt"))
        save_cnv_ids(rare_10kb[CNV_ID_COL], os.path.join(output_dir, f"{anc}_rare_10kb_cnvs.txt"))
        save_cnv_ids(genic[CNV_ID_COL],     os.path.join(output_dir, f"{anc}_genic_rare_cnvs.txt"))

        print(
            f"[{anc}]  All rare: {rare[CNV_ID_COL].nunique():>7,} | "
            f">=10kb: {rare_10kb[CNV_ID_COL].nunique():>6,} | "
            f"Genic: {genic[CNV_ID_COL].nunique():>7,}"
        )


# =============================================================================
# MODULE 4: PLOT CNV COUNT SUMMARIES
# =============================================================================

def plot_cnv_counts(annotations: dict, plot_dir: str) -> None:
    """Bar charts: All rare / Genic / >=10kb / <10kb per ancestry."""
    os.makedirs(plot_dir, exist_ok=True)
    for anc, df in annotations.items():
        rare      = filter_rare_del_dup(df)
        rare_10kb = filter_rare_del_dup(df, min_size_bp=SIZE_10KB)

        rare[OVERLAP_COL] = pd.to_numeric(rare.get(OVERLAP_COL, pd.Series(dtype=float)), errors="coerce")
        genic_mask = (
            (rare[GENE_COL].notna() & (rare[GENE_COL] != "."))
            | (rare[OVERLAP_COL] >= 1)
        )

        counts = {
            "All Rare":      rare[CNV_ID_COL].nunique(),
            "Genic":         rare[genic_mask][CNV_ID_COL].nunique(),
            "Large (>=10kb)": rare_10kb[CNV_ID_COL].nunique(),
            "Small (<10kb)": rare[rare[LENGTH_COL] < SIZE_10KB][CNV_ID_COL].nunique(),
        }

        fig, ax = plt.subplots(figsize=(8, 6))
        bars = ax.bar(counts.keys(), counts.values(),
                      color=["skyblue", "salmon", "lightgreen", "orange"])

        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2,
                    h + max(counts.values()) * 0.01,
                    f"{int(h):,}", ha="center", va="bottom", fontweight="bold")

        ax.set_title(f"{anc}: Unique Rare CNV Counts")
        ax.set_ylabel("Number of Unique IDs")
        ax.tick_params(axis="x", rotation=15)
        ax.grid(axis="y", linestyle="--", alpha=0.7)
        fig.text(0.5, -0.05,
                 "Note: Genic = gene overlap ≥1bp. Autosomes only.",
                 ha="center", fontsize=9, style="italic")

        out = os.path.join(plot_dir, f"{anc}_cnv_count_summary.png")
        fig.savefig(out, bbox_inches="tight", dpi=300)
        plt.close(fig)
        print(f"[{anc}] Plot saved: {out}")


# =============================================================================
# MODULE 5: GENE-SET OVERLAP EXTRACTION (UNIVERSAL)
# =============================================================================

def load_gene_sets_from_multi_set_file(path: str) -> dict:
    """
    Load a file with 'Set' and 'gene' columns.
    Returns {set_name: frozenset(genes)}.
    """
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lower().strip() for c in df.columns]
    df = df.dropna(subset=["gene", "set"])
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    return {
        name: frozenset(grp["gene"])
        for name, grp in df.groupby("set")
    }


def load_gene_sets_from_file_list(file_paths: list) -> dict:
    """
    Each file is a single-column gene list (no header).
    Returns {filename_stem: frozenset(genes)}.
    """
    result = {}
    for path in file_paths:
        if not os.path.exists(path):
            print(f"  WARNING: not found — {path}")
            continue
        stem = os.path.basename(path).replace(".txt", "")
        df = pd.read_csv(path, sep="\t", header=None, names=["gene"])
        result[stem] = frozenset(
            df["gene"].dropna().astype(str).str.upper().str.strip()
        )
    return result


def load_gene_sets_from_directory(directory: str) -> dict:
    """
    Load every *.txt file in a directory as a gene list.
    Returns {filename_stem: frozenset(genes)}.
    """
    result = {}
    for fname in sorted(os.listdir(directory)):
        if not fname.endswith(".txt"):
            continue
        stem = fname.replace(".txt", "")
        df = pd.read_csv(
            os.path.join(directory, fname),
            sep="\t", header=None, names=["gene"]
        )
        result[stem] = frozenset(
            df["gene"].dropna().astype(str).str.upper().str.strip()
        )
    return result


def extract_cnv_ids_for_gene_sets(
    annotations: dict,
    gene_sets: dict,
    output_subdir: str,
    min_size_bp: int = 0,
    file_prefix: str = "cnv_ids",
) -> pd.DataFrame:
    """
    For each (ancestry, gene_set) combination, find overlapping rare CNV IDs
    and write them to individual files.  Returns a combined summary DataFrame.

    Parameters
    ----------
    annotations   : {ancestry: DataFrame}
    gene_sets     : {set_name: frozenset(genes)}
    output_subdir : base output directory (ancestry subfolder created automatically)
    min_size_bp   : 0 = all sizes; 10_000 = >=10kb only
    file_prefix   : prefix for output filenames
    """
    size_label = f"_{min_size_bp // 1000}kb" if min_size_bp else ""
    all_summary = []

    for anc, df in annotations.items():
        out_dir = os.path.join(output_subdir, anc)
        os.makedirs(out_dir, exist_ok=True)

        rare      = filter_rare_del_dup(df, min_size_bp=min_size_bp)
        cnv_genes = expand_genes(rare)

        print(f"\n[{anc}] Rare{size_label} DEL/DUP: {rare[CNV_ID_COL].nunique():,} | "
              f"Genic: {cnv_genes[CNV_ID_COL].nunique():,}")

        for set_name, genes in gene_sets.items():
            cnv_ids, overlap_genes = get_overlapping(cnv_genes, genes)

            fname = f"{file_prefix}{size_label}_{safe_name(set_name)}.txt"
            save_cnv_ids(cnv_ids, os.path.join(out_dir, fname))

            all_summary.append({
                "Ancestry":          anc,
                "Set":               set_name,
                "Min_size_bp":       min_size_bp,
                "Total_Genes_Input": len(genes),
                "Overlapping_Genes": overlap_genes.nunique(),
                "Total_CNVs":        cnv_ids.nunique(),
            })

            print(f"  {set_name:<50}  CNVs: {cnv_ids.nunique():>6,}  "
                  f"genes: {overlap_genes.nunique()}/{len(genes)}")

        # Save per-ancestry summary
        summary = pd.DataFrame([r for r in all_summary if r["Ancestry"] == anc])
        summary_path = os.path.join(
            out_dir,
            f"{anc}{size_label}_{os.path.basename(output_subdir)}_summary.tsv"
        )
        summary.to_csv(summary_path, sep="\t", index=False)

    return pd.DataFrame(all_summary)


# =============================================================================
# MODULE 6: LENGTH MAP (UTILITY)
# =============================================================================

def build_cnv_length_map(
    annotation_file: str,
    output_file: str = "cnv_length_map.tsv",
    id_col: str = "ID",
) -> pd.DataFrame:
    """Extract CNV_ID → total_cnv_length_bp map from a sites-only annotation file."""
    df = pd.read_csv(annotation_file, sep="\t", usecols=[id_col, LENGTH_COL])
    df = df.rename(columns={id_col: "CNV_ID"}).drop_duplicates(subset=["CNV_ID"])
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Length map saved: {output_file}  ({len(df):,} unique CNV IDs)")
    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CNV BURDEN ANALYSIS PIPELINE")
    print("=" * 70)

    # -- Load annotations ---------------------------------------------------
    print("\n[1/6] Loading annotations...")
    annotations = load_annotations(ANNOTATION_FILES)
    for anc, df in annotations.items():
        print(f"  {anc}: {len(df):,} rows")

    # -- Export rare annotations --------------------------------------------
    print("\n[2/6] Exporting rare annotations...")
    export_rare_annotations(annotations, "./cnv_vcf_plink/annotation")

    # -- Size statistics ----------------------------------------------------
    print("\n[3/6] Size statistics...")
    stats = print_size_stats(annotations)

    # -- Export CNV ID lists ------------------------------------------------
    print("\n[4/6] Exporting CNV ID lists (all / 10kb / genic)...")
    export_cnv_id_lists(annotations, "./cnv_vcf_plink")

    # -- Plots --------------------------------------------------------------
    print("\n[5/6] Generating plots...")
    plot_cnv_counts(annotations, PLOT_DIR)

    # -- Gene-set overlaps --------------------------------------------------
    print("\n[6/6] Gene-set overlap extraction...")

    # AUD gene sets (multi-set file with Set + gene columns)
    aud_sets = load_gene_sets_from_multi_set_file(GENE_SET_FILES["aud"])
    print(f"  AUD gene sets loaded: {len(aud_sets)}")

    for min_bp, label in [(0, "aud"), (SIZE_10KB, "aud_10kb")]:
        extract_cnv_ids_for_gene_sets(
            annotations, aud_sets,
            output_subdir=os.path.join(OUTPUT_ROOT, label),
            min_size_bp=min_bp,
            file_prefix="aud_cnv_ids",
        )

    # NDD gene sets (multi-set file with Set + gene columns)
    ndd_sets = load_gene_sets_from_multi_set_file(GENE_SET_FILES["ndd"])
    print(f"  NDD gene sets loaded: {len(ndd_sets)}")

    for min_bp, label in [(0, "ndd"), (SIZE_10KB, "ndd_10kb")]:
        extract_cnv_ids_for_gene_sets(
            annotations, ndd_sets,
            output_subdir=os.path.join(OUTPUT_ROOT, label),
            min_size_bp=min_bp,
            file_prefix="ndd_pgc_cnv_ids",
        )

    # PGC GWAS gene sets (individual files, no header)
    pgc_sets = load_gene_sets_from_file_list(GENE_SET_FILES["pgc_gwas"])
    print(f"  PGC GWAS gene sets loaded: {len(pgc_sets)}")

    for min_bp, label in [(0, "pgc_gwas"), (SIZE_10KB, "pgc_gwas_10kb")]:
        extract_cnv_ids_for_gene_sets(
            annotations, pgc_sets,
            output_subdir=os.path.join(OUTPUT_ROOT, label),
            min_size_bp=min_bp,
            file_prefix="pgc_gwas_cnv_ids",
        )

    # Abnormal gene sets (directory of single-column files)
    abnormal_sets = load_gene_sets_from_directory(GENE_SET_FILES["abnormal_dir"])
    print(f"  Abnormal gene sets loaded: {len(abnormal_sets)}")

    for min_bp, label in [(0, "abnormal"), (SIZE_10KB, "abnormal_10kb")]:
        extract_cnv_ids_for_gene_sets(
            annotations, abnormal_sets,
            output_subdir=os.path.join(OUTPUT_ROOT, label),
            min_size_bp=min_bp,
            file_prefix="abnormal_cnv_ids",
        )

    # Housekeeping / negative control (multi-set file with Set + gene columns)
    hk_sets = load_gene_sets_from_multi_set_file(GENE_SET_FILES["housekeeping"])
    print(f"  Housekeeping gene sets loaded: {len(hk_sets)}")

    for min_bp, label in [(0, "housekeeping"), (SIZE_10KB, "housekeeping_10kb")]:
        extract_cnv_ids_for_gene_sets(
            annotations, hk_sets,
            output_subdir=os.path.join(OUTPUT_ROOT, label),
            min_size_bp=min_bp,
            file_prefix="housekeeping_cnv_ids",
        )

    print("\n✓ Pipeline complete.")


if __name__ == "__main__":
    main()