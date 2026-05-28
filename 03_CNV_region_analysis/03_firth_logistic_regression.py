import os 
import pandas as pd
import subprocess
import time

my_bucket = os.getenv('WORKSPACE_BUCKET')

#####################################################################################
# Step 3: Run Firth's Logistic Regression
# This script performs the association analysis for all ancestry groups and CNV types by:
# 1. Syncing covariate tables to the Google Cloud Bucket.
# 2. Submitting dsub jobs to run PLINK 2.0 --glm firth across all chromosomes.
#####################################################################################

# Define the analysis groups (Ancestry, Analysis Type) and copy the covariate tables
groups = [
    ('EUR', 'WITHIN'),
    ('AFR', 'WITHIN'),
    ('AMR', 'WITHIN'),
    ('EAS', 'WITHIN'),
    ('CSA', 'WITHIN'),
    ('MID', 'WITHIN'),
    ('TRANS', 'TRANS')
]

# Paths
plink_data_dir = f"{my_bucket}/data/cnv_vcf_plink/plink_files"
local_covar_dir = "./aud/AUD_Cases_Controls_EHR_Depth_Within_Ancestry"

# 1. Sync covariate tables from local directory to the cloud bucket
for anc, analysis in groups:
    covar_file = f"all_covariates_{analysis}_{anc}_AUD_age_sex_pcs.txt"
    cloud_dest = f"{my_bucket}/data/aud/AUD_Cases_Controls_EHR_Depth_{analysis}_Ancestry_analysis/"
    subprocess.run(["gsutil", "cp", f"{local_covar_dir}/{covar_file}", cloud_dest], check=True)

# 2. Gather list of PLINK file prefixes (using .bed files to identify chromosomes)
ls_cmd = f"gsutil -u $GOOGLE_PROJECT ls {plink_data_dir}/*.bed" 
file_list_gs = subprocess.check_output(ls_cmd, shell=True)
decoded_file_list =  file_list_gs.decode('utf-8')
file_list_dir = decoded_file_list.strip().split('\n')
prefixes = [item.split('/')[-1].replace(".bed", "") for item in file_list_dir]

# 3. Submit Jobs for each Ancestry and CNV type
for anc, analysis in groups:
    covar_dir = f"{my_bucket}/data/aud/AUD_Cases_Controls_EHR_Depth_{analysis}_Ancestry_analysis"
    covar_file = f"all_covariates_{analysis}_{anc}_AUD_age_sex_pcs.txt"
    anc_dir = f"{my_bucket}/data/aud/AUD_Cases_Controls_EHR_Depth_{analysis}_Ancestry_analysis/{analysis}_{anc}"

    for prefix in prefixes:
        for cnv in ['DELs', 'DUPs']:
            job_name = f"firth_{anc}_{prefix}_{cnv}"
            
            dsub_command = f'''
            source ~/aou_dsub.bash  

            aou_dsub \\
                --image biocontainer/plink2:alpha2.3_jan2020 \\
                --name "{job_name}" \\
                --boot-disk-size 100 \\
                --disk-size 100 \\
                --logging "{anc_dir}/logs/{cnv}/" \\
                --input-recursive plink_input="{plink_data_dir}" \\
                --input-recursive cov_input="{covar_dir}" \\
                --output-recursive output="{anc_dir}/results/{cnv}/" \\
                --command 'plink2 \\
                            --bfile ${{plink_input}}/{prefix} \\
                            --glm firth \\
                            --out ${{output}}/{prefix}_{anc}_{cnv}_firth \\
                            --covar ${{cov_input}}/{covar_file} \\
                            --covar-variance-standardize \\
                            --extract ${{plink_input}}/{prefix}.bim_{cnv}.txt' 
            '''
            # Submitting the job to Google Cloud Batch
            subprocess.run(dsub_command, shell=True, executable='/bin/bash')
            print(f'Job submitted: {job_name}')
            time.sleep(0.5)

print("Association testing submission complete.")