"""
================================================================================
  UNIVERSAL CNV BURDEN ANALYSIS PIPELINE
================================================================================
  *** NOTE: PLEASE READ THIS ENTIRE HEADER BEFORE RUNNING THE SCRIPT.  ***

  Filters, annotates, and extracts rare CNV IDs overlapping user-defined gene
  sets for downstream burden analyses (e.g. PLINK --extract).

  USAGE
  -----
  Edit the CONFIG block below, then run:
      python cnv_burden_pipeline.py

  INPUT REQUIREMENTS
  ------------------
  1. Master CNV annotation TSV  (--master_annotation)
       Required columns: ID, MAF, chr (or chr_2), total_cnv_length_bp,
                         gene, cnv_gene_overlap_length
  2. Gene set files (--gene_sets_dir  OR  --gene_set_files)
       Accepted layouts:
         • Single-column, no header  →  one gene per line
         • Two-column TSV with header  →  columns "gene" + "Set"
             (multiple sets in one file are split automatically)
  3. (Optional) Pre-built rare-10 kb master list  (--rare_10kb_ids_file)
       Plain text, one CNV ID per line.  If not supplied the script derives
       it automatically from the annotation file.

  OUTPUT
  ------
  All files are written to --output_dir (default: ./pipeline_output).

  Per-gene-set files
  ~~~~~~~~~~~~~~~~~~
    {set_name}_all_rare.txt          IDs of all rare (MAF < threshold) DEL/DUP
    {set_name}_rare_10kb.txt         IDs of rare + ≥ min_size_bp DEL/DUP

  Shared helper files
  ~~~~~~~~~~~~~~~~~~~
    rare_del_dup_sites.tsv           Filtered site list (MAF < threshold, DEL/DUP)
    cnv_length_map.tsv               CNV_ID → total_cnv_length_bp lookup
    summary_table.csv                Per-set counts (all rare / rare-large)

  DEPENDENCIES
  ------------
      pip install pandas
================================================================================
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# INLINE CONFIG  ←  edit this block instead of using command-line arguments
# (command-line arguments override these values when supplied)
# ──────────────────────────────────────────────────────────────────────────────
CONFIG = dict(
    # ── Required inputs ───────────────────────────────────────────────────────
    master_annotation  = "AoU_srWGS_SV.v8.sites_only.vcf_Extracted_ID_POS_SVLEN_AF_cytoband_genes.tsv",

    # Gene sets: supply a directory, a list of files, or a two-column TSV.
    # All three options can be set simultaneously; results are merged.
    gene_sets_dir      = "GENE-SETS",          # all *.txt files in this folder
    gene_set_files     = [                      # individual files
        # "GENE-SETS/pgc_ptsd_cnv_ndd_gene_sets.txt",
        # "GENE-SETS/gene_housekeeping_set_neg_control.txt",
    ],

    # ── Optional pre-built filter list ────────────────────────────────────────
    # Leave as "" to auto-derive from the annotation file.
    rare_10kb_ids_file = "",

    # ── Filtering parameters ──────────────────────────────────────────────────
    maf_threshold      = 0.01,   # keep CNVs with MAF < this value
    min_size_bp        = 10_000, # "large" CNV threshold (bp)
    sv_types           = ["DEL", "DUP"],  # SV types to include (matched against ID)

    # ── Sex-chromosome handling ────────────────────────────────────────────────
    # Columns that may hold chromosome information
    chr_columns        = ["chr", "chr_2"],
    sex_chromosomes    = ["chrX", "chrY", "X", "Y"],
    autosomes_only     = False,   # set True to drop sex-chromosome CNVs

    # ── Output ────────────────────────────────────────────────────────────────
    output_dir         = "pipeline_output",
    log_level          = "INFO",   # DEBUG | INFO | WARNING | ERROR
)
# ──────────────────────────────────────────────────────────────────────────────


# ══════════════════════════════════════════════════════════════════════════════
# LOGGING
# ══════════════════════════════════════════════════════════════════════════════

def _setup_logging(level: str) -> logging.Logger:
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
        level=getattr(logging, level.upper(), logging.INFO),
    )
    return logging.getLogger("cnv_pipeline")


# ══════════════════════════════════════════════════════════════════════════════
# GENE-SET LOADING
# ══════════════════════════════════════════════════════════════════════════════

def _infer_layout(path: Path) -> str:
    """Return 'multi' if file has header with 'gene' + 'Set', else 'single'."""
    try:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.str.lower().tolist()
        if "gene" in header and "set" in header:
            return "multi"
    except Exception:
        pass
    return "single"


def load_gene_sets(
    files: list[Path],
    log: logging.Logger,
) -> dict[str, set[str]]:
    """
    Load an arbitrary mix of single-column and two-column gene-set files.

    Returns
    -------
    dict mapping set_name → frozenset of upper-cased gene symbols.
    """
    sets: dict[str, set[str]] = {}

    for path in files:
        if not path.exists():
            log.warning("Gene-set file not found, skipping: %s", path)
            continue

        layout = _infer_layout(path)

        if layout == "multi":
            df = pd.read_csv(path, sep="\t")
            df.columns = df.columns.str.lower().str.strip()
            df = df[["gene", "set"]].dropna()
            df["gene"] = df["gene"].str.upper().str.strip()
            for set_name, grp in df.groupby("set"):
                genes = set(grp["gene"])
                sets.setdefault(str(set_name), set()).update(genes)
            log.info("  %s  →  multi-set layout, found sets: %s",
                     path.name, sorted(sets.keys()))
        else:
            
            set_name = path.stem
            df = pd.read_csv(path, sep="\t", header=None, names=["gene"])
            df = df["gene"].dropna().astype(str).str.upper().str.strip()
            genes = set(df.unique())
            sets.setdefault(set_name, set()).update(genes)
            log.info("  %s  →  single-column layout, %d genes", path.name, len(genes))

    total = sum(len(v) for v in sets.values())
    log.info("Gene sets loaded: %d sets, %d gene entries total", len(sets), total)
    return sets


def _collect_gene_set_paths(cfg: dict) -> list[Path]:
    """Gather every gene-set file from the config."""
    paths: list[Path] = []
   
    for f in cfg.get("gene_set_files") or []:
        p = Path(f)
        if p.is_file():
            paths.append(p)

    d = cfg.get("gene_sets_dir", "")
    if d and Path(d).is_dir():
        found = sorted(Path(d).glob("*.txt"))
        paths.extend(found)

    return list(dict.fromkeys(paths)) 


# ══════════════════════════════════════════════════════════════════════════════
# CNV ANNOTATION LOADING & FILTERING
# ══════════════════════════════════════════════════════════════════════════════

def load_master_annotation(path: str, log: logging.Logger) -> pd.DataFrame:
    log.info("Loading master annotation: %s", path)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    log.info("  Rows loaded: %d", len(df))
    return df


def filter_rare_cnvs(
    df: pd.DataFrame,
    maf_threshold: float,
    sv_types: list[str],
    chr_columns: list[str],
    sex_chromosomes: list[str],
    autosomes_only: bool,
    log: logging.Logger,
) -> pd.DataFrame:
    """
    Apply MAF, SV-type, and optional autosome filters.
    Returns the filtered DataFrame.
    """
    # ── MAF filter ────────────────────────────────────────────────────────────
    df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")
    before = len(df)
    df = df[df["MAF"] < maf_threshold].copy()
    log.info("After MAF < %.3f filter: %d → %d rows", maf_threshold, before, len(df))

    # ── SV type filter (match against ID column) ──────────────────────────────
    pattern = "|".join(sv_types)
    before = len(df)
    df = df[df["ID"].astype(str).str.contains(pattern, na=False)].copy()
    log.info(
        "After SV-type filter (%s): %d → %d rows",
        "/".join(sv_types), before, len(df),
    )

    # ── Autosome filter (optional) ────────────────────────────────────────────
    if autosomes_only:
        chr_col = next((c for c in chr_columns if c in df.columns), None)
        if chr_col:
            before = len(df)
            df = df[~df[chr_col].astype(str).str.strip().isin(sex_chromosomes)].copy()
            log.info(
                "After autosome-only filter (via '%s'): %d → %d rows",
                chr_col, before, len(df),
            )
        else:
            log.warning(
                "autosomes_only=True but none of %s found in columns; skipping.",
                chr_columns,
            )

    log.info("Rare DEL/DUP sites: %d unique IDs", df["ID"].nunique())
    return df


# ══════════════════════════════════════════════════════════════════════════════
# LARGE-CNV MASTER LIST
# ══════════════════════════════════════════════════════════════════════════════

def get_large_cnv_ids(
    rare_df: pd.DataFrame,
    min_size_bp: int,
    ids_file: str,
    log: logging.Logger,
) -> set[str]:
    """
    Return the set of rare CNV IDs that are ≥ min_size_bp.
    If ids_file is supplied and exists, load from there instead.
    """
    if ids_file and Path(ids_file).exists():
        log.info("Loading pre-built large-CNV list: %s", ids_file)
        with open(ids_file) as fh:
            ids = {line.strip() for line in fh if line.strip()}
        log.info("  Loaded %d IDs from file.", len(ids))
        return ids

    log.info("Deriving large-CNV list (>= %d bp) from annotation …", min_size_bp)
    large = rare_df[rare_df["total_cnv_length_bp"] >= min_size_bp]["ID"]
    ids = set(large.unique())
    log.info("  Derived %d large CNV IDs.", len(ids))
    return ids


# ══════════════════════════════════════════════════════════════════════════════
# GENIC OVERLAP
# ══════════════════════════════════════════════════════════════════════════════

def build_gene_id_index(rare_df: pd.DataFrame, log: logging.Logger) -> pd.DataFrame:
    """
    Return a slim DataFrame with columns [gene_upper, ID] representing
    which CNV IDs overlap which genes.
    """
    required = {"gene", "ID", "cnv_gene_overlap_length"}
    missing = required - set(rare_df.columns)
    if missing:
        log.warning(
            "Columns %s missing; genic filter may be incomplete.", missing
        )

    genic_mask = pd.Series(False, index=rare_df.index)
    if "gene" in rare_df.columns:
        genic_mask |= rare_df["gene"].astype(str).str.strip() != "."
    if "cnv_gene_overlap_length" in rare_df.columns:
        genic_mask |= pd.to_numeric(
            rare_df["cnv_gene_overlap_length"], errors="coerce"
        ).fillna(0) >= 1

    genic_df = rare_df[genic_mask][["gene", "ID"]].copy()
    genic_df["gene_upper"] = genic_df["gene"].astype(str).str.upper().str.strip()
    log.info(
        "Genic subset: %d rows, %d unique CNV IDs",
        len(genic_df), genic_df["ID"].nunique(),
    )
    return genic_df[["gene_upper", "ID"]]


# ══════════════════════════════════════════════════════════════════════════════
# PER-SET EXTRACTION
# ══════════════════════════════════════════════════════════════════════════════

def extract_ids_for_set(
    set_name: str,
    set_genes: set[str],
    gene_id_index: pd.DataFrame,
    large_ids: set[str],
    out_dir: Path,
    log: logging.Logger,
) -> dict[str, int]:
    """
    Write two ID lists for a single gene set.
    Returns counts for the summary table.
    """
    safe_name = set_name.replace(" ", "_").replace("/", "-")

    # All rare genic CNVs
    all_rare_ids = set(
        gene_id_index[gene_id_index["gene_upper"].isin(set_genes)]["ID"].unique()
    )

    # Large (≥ min_size_bp) subset
    rare_large_ids = all_rare_ids & large_ids

    def write_ids(ids: set, fname: str) -> None:
        (out_dir / fname).write_text("\n".join(sorted(ids)))

    write_ids(all_rare_ids,   f"{safe_name}_all_rare.txt")
    write_ids(rare_large_ids, f"{safe_name}_rare_large.txt")

    log.info(
        "%-50s  all_rare=%d  rare_large=%d",
        set_name, len(all_rare_ids), len(rare_large_ids),
    )
    return {
        "Gene_Set":        set_name,
        "Genes_in_Set":    len(set_genes),
        "All_Rare_CNVs":   len(all_rare_ids),
        "Rare_Large_CNVs": len(rare_large_ids),
    }


# ══════════════════════════════════════════════════════════════════════════════
# SHARED HELPER FILES
# ══════════════════════════════════════════════════════════════════════════════

def write_shared_files(
    rare_df: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger,
) -> None:
    """Write the filtered site list and the ID→length lookup."""

    # Filtered site list
    site_path = out_dir / "rare_del_dup_sites.tsv"
    rare_df.to_csv(site_path, sep="\t", index=False)
    log.info("Filtered site list → %s", site_path)

    # Length map
    if "total_cnv_length_bp" in rare_df.columns:
        length_map = (
            rare_df[["ID", "total_cnv_length_bp"]]
            .drop_duplicates(subset=["ID"])
            .rename(columns={"ID": "CNV_ID"})
        )
        lmap_path = out_dir / "cnv_length_map.tsv"
        length_map.to_csv(lmap_path, sep="\t", index=False)
        log.info("Length map → %s  (%d unique IDs)", lmap_path, len(length_map))


# ══════════════════════════════════════════════════════════════════════════════
# REPORTING
# ══════════════════════════════════════════════════════════════════════════════

def print_and_save_summary(
    summary_rows: list[dict],
    out_dir: Path,
    log: logging.Logger,
) -> None:
    df = (
        pd.DataFrame(summary_rows)
        .sort_values("All_Rare_CNVs", ascending=False)
        .reset_index(drop=True)
    )
    csv_path = out_dir / "summary_table.csv"
    df.to_csv(csv_path, index=False)
    log.info("Summary table → %s", csv_path)

    # Console table
    col_w = max(len(r["Gene_Set"]) for r in summary_rows) + 2
    header = (
        f"\n{'GENE SET':<{col_w}} | {'GENES':>8} | "
        f"{'ALL RARE':>10} | {'RARE LARGE':>12}"
    )
    sep = "─" * len(header)
    print(sep)
    print(header)
    print(sep)
    for _, row in df.iterrows():
        print(
            f"{row['Gene_Set']:<{col_w}} | "
            f"{row['Genes_in_Set']:>8,} | "
            f"{row['All_Rare_CNVs']:>10,} | "
            f"{row['Rare_Large_CNVs']:>12,}"
        )
    print(sep)


# ══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSER  (values override CONFIG)
# ══════════════════════════════════════════════════════════════════════════════

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Universal CNV Burden Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--master_annotation",  default=CONFIG["master_annotation"])
    p.add_argument("--gene_sets_dir",      default=CONFIG["gene_sets_dir"])
    p.add_argument("--gene_set_files",     nargs="*", default=CONFIG["gene_set_files"])
    p.add_argument("--rare_10kb_ids_file", default=CONFIG["rare_10kb_ids_file"])
    p.add_argument("--maf_threshold",      type=float, default=CONFIG["maf_threshold"])
    p.add_argument("--min_size_bp",        type=int,   default=CONFIG["min_size_bp"])
    p.add_argument("--sv_types",           nargs="+",  default=CONFIG["sv_types"])
    p.add_argument("--autosomes_only",     action="store_true",
                   default=CONFIG["autosomes_only"])
    p.add_argument("--output_dir",         default=CONFIG["output_dir"])
    p.add_argument("--log_level",          default=CONFIG["log_level"])
    return p


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    log  = _setup_logging(args.log_level)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", out_dir.resolve())

    # ── 1. Load master annotation ─────────────────────────────────────────────
    master_df = load_master_annotation(args.master_annotation, log)

    # ── 2. Filter rare DEL/DUP CNVs ──────────────────────────────────────────
    rare_df = filter_rare_cnvs(
        master_df,
        maf_threshold  = args.maf_threshold,
        sv_types       = args.sv_types,
        chr_columns    = CONFIG["chr_columns"],
        sex_chromosomes= CONFIG["sex_chromosomes"],
        autosomes_only = args.autosomes_only,
        log            = log,
    )

    # ── 3. Write shared helper files ──────────────────────────────────────────
    write_shared_files(rare_df, out_dir, log)

    # ── 4. Build large-CNV master set ─────────────────────────────────────────
    large_ids = get_large_cnv_ids(
        rare_df,
        min_size_bp = args.min_size_bp,
        ids_file    = args.rare_10kb_ids_file,
        log         = log,
    )

    # ── 5. Build gene → ID index ──────────────────────────────────────────────
    gene_id_index = build_gene_id_index(rare_df, log)

    # ── 6. Load all gene sets ─────────────────────────────────────────────────
    gs_paths = _collect_gene_set_paths(vars(args))
    if not gs_paths:
        log.error(
            "No gene-set files found. "
            "Set gene_sets_dir or gene_set_files in CONFIG or pass --gene_sets_dir."
        )
        sys.exit(1)

    gene_sets = load_gene_sets(gs_paths, log)
    if not gene_sets:
        log.error("Gene sets dictionary is empty after loading. Check file formats.")
        sys.exit(1)

    # ── 7. Per-set extraction ─────────────────────────────────────────────────
    log.info("=" * 60)
    log.info("Extracting CNV IDs per gene set …")
    log.info("=" * 60)

    summary_rows: list[dict] = []
    for set_name, set_genes in sorted(gene_sets.items()):
        row = extract_ids_for_set(
            set_name      = set_name,
            set_genes     = set_genes,
            gene_id_index = gene_id_index,
            large_ids     = large_ids,
            out_dir       = out_dir,
            log           = log,
        )
        summary_rows.append(row)

    # ── 8. Summary ────────────────────────────────────────────────────────────
    if summary_rows:
        print_and_save_summary(summary_rows, out_dir, log)

    log.info("Pipeline complete.  All outputs in: %s", out_dir.resolve())


if __name__ == "__main__":
    main()