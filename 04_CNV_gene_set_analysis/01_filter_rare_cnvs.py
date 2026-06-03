"""
filter_rare_cnvs.py
-------------------
Identifies and extracts rare/common CNVs (deletions and duplications) from
PLINK .bim files using an allele frequency annotation BED file.

Usage:
    python filter_rare_cnvs.py [--bim-dir DIR] [--af-file FILE] \
                               [--n-samples N] [--maf-max FLOAT] \
                               [--mac-thresholds INT [INT ...]] \
                               [--out-dir DIR]

Defaults match the original hardcoded paths.
"""

import argparse
import os
import warnings

import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bim-dir",  default="./cnv_vcf_plink/SV_bim_files/",
                   help="Directory containing PLINK .bim files.")
    p.add_argument("--af-file",
                   default="./cnv_vcf_plink/AoU_srWGS_SV.v8_sites_only_coordinates_AF_cytoband_unique.bed",
                   help="7-column BED file with CNV allele frequencies.")
    p.add_argument("--out-dir",  default=None,
                   help="Output directory (defaults to --bim-dir).")
    p.add_argument("--n-samples", type=int, default=97_061,
                   help="Total number of samples used to compute AF.")
    p.add_argument("--maf-max",  type=float, default=0.01,
                   help="MAF ceiling for 'rare' CNVs (default: 0.01).")
    p.add_argument("--mac-thresholds", type=int, nargs="+", default=[21, 11],
                   help="One or more MAC floor values to run (default: 21 11).")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_af_data(af_file: str) -> pd.DataFrame:
    """Load and clean the CNV allele-frequency annotation BED file."""
    df = pd.read_csv(
        af_file, sep="\t", header=None,
        names=["chrom", "start", "end", "ID", "CNV_type", "AF", "cytoband"],
    )
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    df = df.dropna(subset=["AF"])
    return df


def load_bim_files(bim_dir: str) -> dict[str, pd.DataFrame]:
    """Load every .bim file in *bim_dir* and return {filename: DataFrame}."""
    bim_data: dict[str, pd.DataFrame] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)   # delim_whitespace
        for fname in sorted(os.listdir(bim_dir)):
            if fname.endswith(".bim"):
                path = os.path.join(bim_dir, fname)
                bim_data[fname] = pd.read_csv(path, delim_whitespace=True, header=None)
    return bim_data


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------

AUTOSOMES = ~(pd.Series(dtype=str).isin(["chrX", "chrY"]))  # used as template


def get_rare_ids(
    af_df: pd.DataFrame,
    af_min: float,
    af_max: float,
) -> tuple[set, set, pd.DataFrame, pd.DataFrame]:
    """
    Return (rare_del_ids, rare_dup_ids, rare_del_df, rare_dup_df) for
    CNVs on autosomes with af_min <= AF < af_max.
    """
    autosome_mask = (~af_df["chrom"].isin(["chrX", "chrY"]) &
                     af_df["AF"].between(af_min, af_max, inclusive="left"))

    rare_dels = af_df[(af_df["CNV_type"] == "<DEL>") & autosome_mask]
    rare_dups = af_df[(af_df["CNV_type"] == "<DUP>") & autosome_mask]
    return set(rare_dels["ID"]), set(rare_dups["ID"]), rare_dels, rare_dups


def get_common_ids(af_df: pd.DataFrame, af_threshold: float) -> tuple[set, set]:
    """Return (common_del_ids, common_dup_ids) for AF >= af_threshold on autosomes."""
    autosome_mask = (~af_df["chrom"].isin(["chrX", "chrY"]) &
                     (af_df["AF"] >= af_threshold))
    common_dels = af_df[(af_df["CNV_type"] == "<DEL>") & autosome_mask]
    common_dups = af_df[(af_df["CNV_type"] == "<DUP>") & autosome_mask]
    return set(common_dels["ID"]), set(common_dups["ID"])


def chrom_label(base_name: str) -> str:
    """Extract a tidy chromosome label from a .bim base name, e.g. 'EUR_chr22_qc' → '22'."""
    parts = base_name.split("chr")
    if len(parts) > 1:
        return parts[-1].split("_")[0].split(".")[0].upper()
    return base_name


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_id_list(ids: set, path: str) -> None:
    with open(path, "w") as fh:
        for cnv_id in sorted(ids):
            fh.write(f"{cnv_id}\n")


def save_bim_subset(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, sep="\t", header=False, index=False)


def print_summary_table(title: str, chrom_counts: dict) -> None:
    cw = (15, 20, 20)
    sep = f"{'─'*cw[0]}{'─'*cw[1]}{'─'*cw[2]}"
    print(f"\n{title}")
    print(sep)
    print(f"{'Chromosome':<{cw[0]}}{'Deletions':<{cw[1]}}{'Duplications':<{cw[2]}}")
    print(sep)
    total_dels = total_dups = 0
    sort_key = lambda x: int(x) if x.isdigit() else (23 if x == "X" else 24 if x == "Y" else 99)
    for chrom in sorted(chrom_counts, key=sort_key):
        d, u = chrom_counts[chrom]["DELs"], chrom_counts[chrom]["DUPs"]
        label = f"chr{chrom}" if chrom.isdigit() or chrom in ("X", "Y") else chrom
        print(f"{label:<{cw[0]}}{d:<{cw[1]}}{u:<{cw[2]}}")
        total_dels += d
        total_dups += u
    print(sep)
    print(f"{'TOTAL':<{cw[0]}}{total_dels:<{cw[1]}}{total_dups:<{cw[2]}}")
    print(sep)


# ---------------------------------------------------------------------------
# Analysis sections
# ---------------------------------------------------------------------------

def report_size_statistics(rare_dels: pd.DataFrame, rare_dups: pd.DataFrame) -> None:
    del_sizes = rare_dels["end"] - rare_dels["start"]
    dup_sizes = rare_dups["end"] - rare_dups["start"]
    all_sizes = pd.concat([del_sizes, dup_sizes])

    print("\n" + "─" * 50)
    print("Size Statistics for Rare CNVs")
    print("─" * 50)
    for label, sizes in [("OVERALL", all_sizes), ("DELETIONS", del_sizes), ("DUPLICATIONS", dup_sizes)]:
        print(f"{label}  (n={len(sizes)})")
        print(f"  Median : {sizes.median() / 1000:.2f} kb")
        print(f"  Mean   : {sizes.mean()   / 1000:.2f} kb")
    print("─" * 50)


def run_mannwhitney(rare_dels: pd.DataFrame, rare_dups: pd.DataFrame) -> None:
    del_sizes = rare_dels["end"] - rare_dels["start"]
    dup_sizes = rare_dups["end"] - rare_dups["start"]
    u_stat, p_val = stats.mannwhitneyu(del_sizes, dup_sizes, alternative="less")

    print("\n" + "─" * 50)
    print("Mann-Whitney U Test  (H₁: median DEL size < median DUP size)")
    print("─" * 50)
    print(f"  Median DEL : {del_sizes.median() / 1000:.2f} kb")
    print(f"  Median DUP : {dup_sizes.median() / 1000:.2f} kb")
    print(f"  U-statistic: {u_stat:.4f}")
    print(f"  P-value    : {p_val:.4e}")
    conclusion = ("Significant — deletions are typically smaller than duplications."
                  if p_val < 0.05 else
                  "Not significant — no evidence deletions are smaller.")
    print(f"  Result     : {conclusion}")
    print("─" * 50)


def filter_and_write_bim_files(
    bim_data: dict,
    del_ids: set,
    dup_ids: set,
    out_dir: str,
    tag: str,
) -> dict:
    """
    For each loaded .bim DataFrame, extract rows matching del_ids / dup_ids,
    write output files with *tag* in the name, and return per-chromosome counts.
    """
    chrom_counts: dict = {}
    for fname, bim_df in bim_data.items():
        base = os.path.splitext(fname)[0]
        chrom = chrom_label(base)

        dels_df = bim_df[bim_df[1].isin(del_ids)]
        dups_df = bim_df[bim_df[1].isin(dup_ids)]

        save_bim_subset(dels_df, os.path.join(out_dir, f"{base}_{tag}_DELs.bim.txt"))
        save_bim_subset(dups_df, os.path.join(out_dir, f"{base}_{tag}_DUPs.bim.txt"))

        chrom_counts[chrom] = {"DELs": len(dels_df), "DUPs": len(dups_df)}
    return chrom_counts


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    out_dir = args.out_dir or args.bim_dir
    os.makedirs(out_dir, exist_ok=True)

    # ── Load data once ─────────────────────────────────────────────────────
    print(f"Loading AF data from: {args.af_file}")
    af_df = load_af_data(args.af_file)

    print(f"Loading .bim files from: {args.bim_dir}")
    bim_data = load_bim_files(args.bim_dir)
    print(f"  Found {len(bim_data)} .bim file(s).\n")

    # ── Simple rare/common split (AF ≤ 0.01) ──────────────────────────────
    simple_del_ids, simple_dup_ids, rare_dels_df, rare_dups_df = get_rare_ids(
        af_df, af_min=0.0, af_max=args.maf_max
    )

    # Master ID lists
    all_rare_ids = simple_del_ids | simple_dup_ids
    save_id_list(all_rare_ids, os.path.join(out_dir, "all_rare_cnvs_to_extract.txt"))
    print(f"Master rare CNV list: {len(all_rare_ids)} IDs saved.")

    save_id_list(simple_del_ids, os.path.join(out_dir, "rare_del_ids_to_extract.txt"))
    save_id_list(simple_dup_ids, os.path.join(out_dir, "rare_dup_ids_to_extract.txt"))

    common_del_ids, common_dup_ids = get_common_ids(af_df, args.maf_max)
    all_common_ids = common_del_ids | common_dup_ids
    save_id_list(all_common_ids, os.path.join(out_dir, "all_common_cnvs_to_extract.txt"))
    save_id_list(common_del_ids, os.path.join(out_dir, "common_del_ids_to_extract.txt"))
    save_id_list(common_dup_ids, os.path.join(out_dir, "common_dup_ids_to_extract.txt"))
    print(f"Master common CNV list: {len(all_common_ids)} IDs saved.")

    # BIM files: simple rare/common split
    filter_and_write_bim_files(bim_data, simple_del_ids, simple_dup_ids, out_dir, "RARE")
    filter_and_write_bim_files(bim_data, common_del_ids, common_dup_ids, out_dir, "COMMON")

    # ── Size statistics & Mann-Whitney (simple rare split) ─────────────────
    report_size_statistics(rare_dels_df, rare_dups_df)
    run_mannwhitney(rare_dels_df, rare_dups_df)

    # ── Dual-threshold filter (MAF < 1% AND MAC ≥ threshold) ───────────────
    for mac_floor in args.mac_thresholds:
        af_min = mac_floor / (2 * args.n_samples)
        tag = f"MAF1_MAC{mac_floor}"
        print(f"\n── Dual filter: MAF < {args.maf_max}  AND  MAC ≥ {mac_floor} "
              f"(AF_min = {af_min:.2e}) ──")

        del_ids, dup_ids, _, _ = get_rare_ids(af_df, af_min=af_min, af_max=args.maf_max)
        print(f"  Rare DELs: {len(del_ids)}   Rare DUPs: {len(dup_ids)}")

        chrom_counts = filter_and_write_bim_files(bim_data, del_ids, dup_ids, out_dir, tag)
        print_summary_table(
            f"Rare CNVs — MAF < {args.maf_max} & MAC ≥ {mac_floor}",
            chrom_counts,
        )

    print("\n✅ All done.")


if __name__ == "__main__":
    main()