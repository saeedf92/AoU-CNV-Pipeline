import pandas as pd
import os

#####################################################################################
# This script updates the 6th column (phenotype) of .fam files for association testing.
# PLINK expects 1 for Control and 2 for Case.
#####################################################################################

####################################################
#### copy .fam files to local workspace (shell) ####
# my_bucket = os.getenv('WORKSPACE_BUCKET')
# ancestries=('EUR' 'AFR' 'AMR' 'EAS' 'CSA' 'MID' 'TRANS')

# # Loop through and copy the files
# for anc in "${ancestries[@]}"; do
#     gsutil -m cp "${my_bucket}/data/cnv_vcf_plink/plink_files/*.fam" "./cohort_within_${anc}_fam_files/"
# done
####################################################

# Define the list of ancestries
ancestries = ['EUR', 'AFR', 'AMR', 'EAS', 'CSA', 'MID', 'TRANS']

# Define the base directory for the files.
base_dir = 'path/to/the/directory'

# Define the list of ancestries
ancestries = ['EUR', 'AFR', 'AMR', 'EAS', 'CSA', 'MID']

# Define the base directory for the files.
base_dir = 'aud/AUD_Cases_Controls_EHR_Depth_Ancestry_PoPMaD/'

def fam_modify(fam_directory, cases_covariate, controls_covariate):
    """
    Modifies the phenotype column of .fam files based on provided cases and controls.

    Args:
        fam_directory (str): The directory containing the .fam files.
        cases_covariate (str): The path to the cases covariate file.
        controls_covariate (str): The path to the controls covariate file.
    """
    # Read AUD cases and controls covariates files
    cases_df = pd.read_csv(cases_covariate, sep="\t")
    controls_df = pd.read_csv(controls_covariate, sep="\t")

    # Ensure IID is treated as strings and trim any whitespace for both cases and controls
    cases_df['IID'] = cases_df['IID'].astype(str).str.strip()
    controls_df['IID'] = controls_df['IID'].astype(str).str.strip()

    # Loop through all .fam files in the directory
    for filename in os.listdir(fam_directory):
        if filename.endswith('.fam'):
            fam_file = os.path.join(fam_directory, filename)
            print(f"Processing {fam_file}...")

            # Read the .fam file
            fam_df = pd.read_csv(fam_file, sep="\t", names=['FID', 'IID', 'IDF', 'IDM', 'SEX', 'pheno'], header=None)

            # Ensure IID is treated as strings and trim any whitespace
            fam_df['IID'] = fam_df['IID'].astype(str).str.strip()

            # First change all pheno to -9
            fam_df['pheno'] = -9

            # Set pheno to 1 for controls
            fam_df.loc[fam_df['IID'].isin(controls_df['IID']), 'pheno'] = 1

            # Set pheno to 2 for cases
            fam_df.loc[fam_df['IID'].isin(cases_df['IID']), 'pheno'] = 2

            # Save the modified .fam file
            fam_df.to_csv(fam_file, sep='\t', header=False, index=False)

    print("All .fam files processed.")

# Loop through each ancestry and call the fam_modify function
for ancestry in ancestries:
    fam_dir = os.path.join(base_dir, f'pheno_{ancestry}_fam_files/')
    cases_file = os.path.join(base_dir, f'covariates_{ancestry}_AUD_cases_age_sex_pcs.txt')
    controls_file = os.path.join(base_dir, f'covariates_{ancestry}_AUD_controls_age_sex_pcs.txt')

    fam_modify(fam_dir, cases_file, controls_file)