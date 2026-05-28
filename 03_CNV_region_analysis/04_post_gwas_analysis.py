import os
import pandas as pd
import numpy as np
import subprocess
from scipy import stats
from statsmodels.stats.multitest import multipletests

my_bucket = os.getenv('WORKSPACE_BUCKET')

#####################################################################################
# Step 4: Post GWAS Analysis
# This script aggregates results, calculates genomic inflation, performs FDR 
# correction, and maps associations to functional annotations.
#####################################################################################

groups = [
    ('EUR', 'WITHIN'),
    ('AFR', 'WITHIN'),
    ('AMR', 'WITHIN'),
    ('EAS', 'WITHIN'),
    ('CSA', 'WITHIN'),
    ('MID', 'WITHIN'),
    ('TRANS', 'TRANS')
]

cnv_types = ['DELs', 'DUPs']

def calculate_lambda(p_values):
    """Calculates the genomic inflation factor lambda."""
    chi2 = stats.chi2.ppf(1 - p_values, 1)
    return np.median(chi2) / stats.chi2.ppf(0.5, 1)

def calculate_lambda_1000(lam, n_cases, n_controls):
    """Calculates lambda standardized to 1000 cases and 1000 controls."""
    return 1 + (lam - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)

def main():
    # Load the Annotation Table generated in Data Processing Step 4
    # Assuming the table is stored locally or in the bucket after Step 4
    annot_path = "aud/annotation_table/CNV_Annotation_Table.csv"
    if os.path.exists(annot_path):
        df_annot = pd.read_csv(annot_path)
        print("Annotation table loaded.")
    else:
        print(f"Warning: Annotation table not found at {annot_path}. Results will not be annotated.")
        df_annot = None

    for anc, analysis in groups:
        for cnv in cnv_types:
            print(f"Processing results for: {anc} - {cnv}")
            
            # Define cloud paths based on Step 3 output structure
            anc_dir = f"{my_bucket}/data/aud/AUD_Cases_Controls_EHR_Depth_{analysis}_Ancestry_analysis/{analysis}_{anc}"
            result_cloud_path = f"{anc_dir}/results/{cnv}/"
            
            # List all .firth files for this group
            ls_cmd = f"gsutil ls {result_cloud_path}*.firth"
            try:
                files = subprocess.check_output(ls_cmd, shell=True).decode('utf-8').strip().split('\n')
            except subprocess.CalledProcessError:
                print(f"  No result files found for {anc} {cnv}. Skipping...")
                continue

            # 4.1 Results Aggregation
            dfs = []
            for f in files:
                # Read results directly from the bucket
                df_temp = pd.read_csv(f, sep='\s+')
                dfs.append(df_temp)
            
            merged_df = pd.concat(dfs, ignore_index=True)
            
            # Filter for valid p-values and the primary additive test
            if 'TEST' in merged_df.columns:
                merged_df = merged_df[merged_df['TEST'] == 'ADD']
            merged_df = merged_df.dropna(subset=['P'])

            # 4.2 Statistical QC
            n_samples = merged_df['OBS_CT'].iloc[0] if 'OBS_CT' in merged_df.columns else 0
            lam = calculate_lambda(merged_df['P'])
            print(f"  Genomic Inflation (Lambda): {lam:.4f}")

            # 4.3 Multiple Testing Correction (FDR)
            _, fdr_p, _, _ = multipletests(merged_df['P'], method='fdr_bh')
            merged_df['FDR_P'] = fdr_p

            # 4.5 Functional Annotation Mapping
            if df_annot is not None:
                # ID column is expected to be consistent across PLINK and Annotation Table
                merged_df = pd.merge(merged_df, df_annot, on='ID', how='left')

            # Save final summary statistics
            output_filename = f"{anc}_{cnv}_summary_stats.csv"
            merged_df.to_csv(output_filename, index=False)
            
            # Upload to bucket
            final_dest = f"{anc_dir}/final_results/"
            subprocess.run(f"gsutil cp {output_filename} {final_dest}", shell=True)
            print(f"  Final summary stats uploaded to {final_dest}")

    print("\nPost-GWAS Analysis Complete.")

if __name__ == "__main__":
    main()