import os
import subprocess
import time
import argparse

# --- GLOBAL CONFIGURATION ---
MY_BUCKET = os.getenv('WORKSPACE_BUCKET')
if not MY_BUCKET:
    raise EnvironmentError("WORKSPACE_BUCKET environment variable not set. Are you running this in the correct environment?")

DATA_DIR = f"{MY_BUCKET}/data/cnv_vcf_plink/plink_files"
COVARIATE_DIR = f"{MY_BUCKET}/data/plink_phenotype_files/plink_CNV_phenotype_01"

# --- CORE JOB SUBMISSION FUNCTION ---
def submit_dsub_job(job_name, log_dir, cnv_list_file, keep_file, output_dir, output_tag, pop):
    """Generates and executes the aou_dsub command for PLINK extraction."""
    dsub_cmd = f"""
    source ~/aou_dsub.bash
    aou_dsub \\
      --image us-central1-docker.pkg.dev/polar-standard-455018-c9/prscs-repo/prscs:1.1.1 \\
      --name "{job_name}" \\
      --boot-disk-size 100 \\
      --disk-size 150 \\
      --machine-type "n2-highmem-8" \\
      --logging "{log_dir}" \\
      --input-recursive plink_input="{DATA_DIR}" \\
      --input cnv_list_file="{cnv_list_file}" \\
      --input keep_file="{keep_file}" \\
      --output-recursive output_dir="{output_dir}" \\
      --command '
        set -euo pipefail
        for bed_file in ${{plink_input}}/AoU_srWGS_SV.v8.chr*.bed; do
            base_prefix=$(basename "${{bed_file}}" .bed)
            raw_output_prefix="${{base_prefix}}{output_tag}"
            echo "Processing ${{base_prefix}} in population {pop}..."
            
            if ! plink2 \\
                --bfile ${{plink_input}}/${{base_prefix}} \\
                --keep ${{keep_file}} \\
                --extract ${{cnv_list_file}} \\
                --recode A \\
                --out ${{output_dir}}/${{raw_output_prefix}}; then
                echo "SKIPPED: No overlapping CNVs found for ${{base_prefix}}"
            fi
        done
      '
    """
    subprocess.run(dsub_cmd, shell=True, executable="/bin/bash")
    time.sleep(0.5)

# --- ANALYSIS MODULES ---

def process_standard(ancestries):
    """Processes standard genic rare CNVs."""
    for pop in ancestries:
        output_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/genic_rare_cnv"
        log_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/logs"
        keep_file = f"{COVARIATE_DIR}/{pop}_standard_code_samplelist.txt"
        cnv_list = f"{MY_BUCKET}/data/gene_sets/{pop}/{pop}_genic_rare_cnvs_to_extract.txt"
        
        print(f"\n[STANDARD] Submitting jobs for {pop}")
        submit_dsub_job(
            job_name=f"GS_{pop}_genic_rCNV_EXTRACT",
            log_dir=log_dir,
            cnv_list_file=cnv_list,
            keep_file=keep_file,
            output_dir=output_dir,
            output_tag=f"_{pop}_genic_rCNV",
            pop=pop
        )

def process_ptsd(ancestries, run_10kb):
    """Processes PTSD GWAS Gene Sets."""
    gene_sets = {
        "pgc_gwas_cnv_ids_protein_coding_suppl_tbl_9.txt": "protein_coding",
        "pgc_gwas_cnv_ids_magma_suppl_tbl_13.txt": "ptsd_magma",
        "pgc_gwas_cnv_ids_brain_tissue_twas_suppl_tbl_16.txt": "brain_tissue_twas", 
        "pgc_gwas_cnv_ids_brain_tissue_eQTL_smr_suppl_tbl_17.txt": "brain_tissue_eQTL_smr",
        "pgc_gwas_cnv_ids_brain_tissue_sc_dlPFC_suppl_tbl_19.txt": "brain_tissue_dlPFC"
    }
    
    folder_type = "pgc_ptsd_gwas_rare_10kb_cnv_count" if run_10kb else "pgc_ptsd_gwas_rare_cnv_count"
    size_tag = "_10kb" if run_10kb else ""

    for pop in ancestries:
        base_out_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/{folder_type}"
        log_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/logs"
        keep_file = f"{COVARIATE_DIR}/{pop}_standard_code_samplelist.txt"

        for file_name, label in gene_sets.items():
            cnv_list = f"{MY_BUCKET}/data/gene_sets/{pop}/ptsd/{file_name}"
            out_dir = f"{base_out_dir}/{file_name.replace('.txt', '')}"
            tag = f"_{pop}_{label}{size_tag}"
            
            print(f"[PTSD] Submitting: {pop} | {label} | 10kb: {run_10kb}")
            submit_dsub_job(f"GS_{tag}", log_dir, cnv_list, keep_file, out_dir, tag, pop)

def process_aud(ancestries, run_10kb):
    """Processes AUD (Alcohol Use Disorder) GWAS Gene Sets."""
    # TODO: Update these keys with your actual AUD gene set file names
    gene_sets = {
        "aud_gwas_cnv_ids_protein_coding_placeholder.txt": "protein_coding",
        "aud_gwas_cnv_ids_magma_placeholder.txt": "aud_magma",
        "aud_gwas_cnv_ids_brain_tissue_placeholder.txt": "brain_tissue"
    }
    
    folder_type = "pgc_aud_gwas_rare_10kb_cnv_count" if run_10kb else "pgc_aud_gwas_rare_cnv_count"
    size_tag = "_10kb" if run_10kb else ""

    for pop in ancestries:
        base_out_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/{folder_type}"
        log_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/logs"
        keep_file = f"{COVARIATE_DIR}/{pop}_standard_code_samplelist.txt"

        for file_name, label in gene_sets.items():
            cnv_list = f"{MY_BUCKET}/data/gene_sets/{pop}/aud/{file_name}"
            out_dir = f"{base_out_dir}/{file_name.replace('.txt', '')}"
            tag = f"_{pop}_{label}{size_tag}"
            
            print(f"[AUD] Submitting: {pop} | {label} | 10kb: {run_10kb}")
            submit_dsub_job(f"GS_{tag}", log_dir, cnv_list, keep_file, out_dir, tag, pop)

def process_ndd(ancestries, run_10kb):
    """Processes Neurodevelopmental Disorder (NDD) Gene Sets."""
    analysis_types = {
        "housekeeping_cnv_ids_present_in_set_HSIAO-HOUSEKEEPING-GENES.txt": "HSIAO_HOUSEKEEPING_GENES",
    }
    
    for pop in ancestries:
        folder_type = "ndd_pgc_rare_cnv_10kb_count" if run_10kb else "ndd_pgc_rare_cnv_count"
        base_out_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/{folder_type}"
        log_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/logs"
        keep_file = f"{COVARIATE_DIR}/{pop}_standard_code_samplelist.txt"

        for base_file, label in analysis_types.items():
            if run_10kb:
                formatted_file = base_file.replace("ndd_pgc_", "ndd_pgc_10kb_").replace("housekeeping_", "housekeeping_10kb_")
                cnv_list = f"{MY_BUCKET}/data/gene_sets/{pop}/ndd/10kb/{formatted_file}"
            else:
                formatted_file = base_file
                cnv_list = f"{MY_BUCKET}/data/gene_sets/{pop}/ndd/{formatted_file}"

            out_dir = f"{base_out_dir}/{formatted_file.replace('.txt', '')}"
            tag = f"_{pop}_{label}"
            
            print(f"[NDD] Submitting: {pop} | {label} | 10kb: {run_10kb}")
            submit_dsub_job(f"GS_{tag}", log_dir, cnv_list, keep_file, out_dir, tag, pop)

def process_abnormal(ancestries, run_10kb, start_idx, end_idx):
    """Processes Abnormal Gene Sets with batching."""
    for pop in ancestries:
        if run_10kb:
            source_dir = f"{MY_BUCKET}/data/gene_sets/{pop}/abnormal/10kb"
            search_pattern = "abnormal_10kb_cnv_ids_present_in_set_*.txt"
            base_out_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/abnormal_gene_sets_rare_10kb_cnv_count"
            replace_str = "abnormal_10kb_cnv_ids_present_in_set_"
        else:
            source_dir = f"{MY_BUCKET}/data/gene_sets/{pop}/abnormal"
            search_pattern = "abnormal_cnv_ids_present_in_set_*.txt"
            base_out_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/rare_cnv/abnormal_gene_sets_rare_cnv_count"
            replace_str = "abnormal_cnv_ids_present_in_set_"

        log_dir = f"{MY_BUCKET}/data/gene_set_output/{pop}/logs"
        keep_file = f"{COVARIATE_DIR}/{pop}_standard_code_samplelist.txt"

        try:
            ls_cmd = f"gsutil ls {source_dir}/{search_pattern}"
            all_files = subprocess.check_output(ls_cmd, shell=True).decode().splitlines()
            all_gene_sets = sorted([os.path.basename(f) for f in all_files])
            batch = all_gene_sets[start_idx:end_idx]
        except subprocess.CalledProcessError:
            print(f"[ERROR] Could not find abnormal sets for {pop}. Skipping.")
            continue

        print(f"\n[ABNORMAL] Processing {pop} | Batch {start_idx}-{end_idx} ({len(batch)} files) | 10kb: {run_10kb}")

        for file_name in batch:
            cnv_list = f"{source_dir}/{file_name}"
            clean_label = file_name.replace(replace_str, "").replace(".txt", "")
            out_dir = f"{base_out_dir}/{clean_label}"
            tag = f"_{pop}_{clean_label}"
            
            submit_dsub_job(f"GS_{pop}_{clean_label}"[:40], log_dir, cnv_list, keep_file, out_dir, tag, pop)

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline to run PLINK extractions for various CNV gene sets.")
    parser.add_argument(
        "--analysis", 
        required=True, 
        choices=["standard", "ptsd", "aud", "ndd", "abnormal"], 
        help="Type of analysis to run."
    )
    parser.add_argument("--ancestries", nargs='+', default=["AMR"], help="List of ancestries (e.g., EUR AFR AMR).")
    parser.add_argument("--run-10kb", action="store_true", help="Flag to use the 10kb filtered sets.")
    parser.add_argument("--start", type=int, default=0, help="Start index for abnormal gene-set batching.")
    parser.add_argument("--end", type=int, default=1000, help="End index for abnormal gene-set batching.")

    args = parser.parse_args()

    print(f"=== Starting CNV Pipeline ===")
    print(f"Target Ancestries: {args.ancestries}")

    if args.analysis == "standard":
        process_standard(args.ancestries)
    elif args.analysis == "ptsd":
        process_ptsd(args.ancestries, args.run_10kb)
    elif args.analysis == "aud":
        process_aud(args.ancestries, args.run_10kb)
    elif args.analysis == "ndd":
        process_ndd(args.ancestries, args.run_10kb)
    elif args.analysis == "abnormal":
        process_abnormal(args.ancestries, args.run_10kb, args.start, args.end)

    print("\n=== Pipeline Submission Complete ===")


#==============================================================================================
#==============================================================================================
# How to use this pipeline EXAMPLE COMMANDS:
# $ python 05_burden_count.py --analysis aud --ancestries EUR AMR --run-10kb

# $ python 05_burden_count.py --analysis ptsd --ancestries EUR AMR --run-10kb

# Run NDD sets for AFR, without 10kb filters:
# $ python 05_burden_count.py --analysis ndd --ancestries AFR

# Run a specific batch of Abnormal sets for AMR:
# $ python 05_burden_count.py --analysis abnormal --ancestries AMR --start 120 --end 147
#==============================================================================================
#==============================================================================================