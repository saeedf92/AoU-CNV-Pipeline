import pandas as pd
import os

my_bucket = os.getenv('WORKSPACE_BUCKET')

#####################################################################################
# Load the cohort file: It should contain FID, IID, and phenotypes columns
#####################################################################################
# The cohort file with FID, IID, and phenotypes
# binary trait: 0/1 for control/case
df_case_control = pd.read_csv('cohort.txt', sep = '\t')

# AoU v8 all participants: it contains FID, IID, age, sex
# age: age at the time of sample collection (continuous variable)
# sex: matching self-reported sex with genetically inferred sex (1 for male, 2 for female, -9 for unknown)
df_aou_v8_all_with_CNV = pd.read_csv('aou_v8_all_participants_with_CNV_data.txt', sep = '\t') 

# This file contains IID, PCs and ancestry information (Singh et al. 2025)
df_ancestry = pd.read_csv('Ancestry.list', sep = '\t')

#####################################################################################
# Merge the cohort file with the ancestry file to get a complete covariate table
#####################################################################################

# Left join age and sex
df_cohort_all_covariates = pd.merge(df_case_control, df_aou_v8_all_with_CNV[['IID', 'age', 'sex']], 
                                    how='left', left_on='IID', right_on='IID')
# Left join PC1-10
### Note: PCs should be added separately based on each ancestry group to create ancestry-specific covariate tables. 
df_cohort_all_covariates = pd.merge(df_cohort_all_covariates, 
                                    df_ancestry[['IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']], 
                                    how='left', left_on='IID', right_on='IID')


# Save the merged dataframe to a new file
df_cohort_all_covariates.to_csv('cohort_covariate_all_covariates_<ancestry_group>.txt', sep = '\t', index=False)
# save cases and controls separately for downstream analyses
df_cohort_all_covariates_cases    = df_cohort_all_covariates[df_cohort_all_covariates['Phenotype'] == 1]
df_cohort_all_covariates_controls = df_cohort_all_covariates[df_cohort_all_covariates['Phenotype'] == 0]
df_cohort_all_covariates_cases.to_csv('cohort_covariate_cases_all_covariates_<ancestry_group>.txt', sep = '\t', index=False)
df_cohort_all_covariates_controls.to_csv('cohort_covariate_controls_all_covariates_<ancestry_group>.txt', sep = '\t', index=False)

#####################################################################################
# Final covariate table should contain the following columns:
# - FID: Family ID (can be the same as IID if not available)
# - IID: Individual ID
# - Phenotype: Binary trait (0 for control, 1 for case)
# - Age: Age at the time of sample collection (continuous variable)
# - Sex: 1 for male, 2 for female, -9 for unknown
# - PCs: Principal components (PC1-10)
#####################################################################################


#####################################################################################################################
# Optimized code for merging and saving covariate tables for cases and controls separately based on ancestry groups.
#####################################################################################################################

# Define the ancestries
ancestries = ['EUR', 'AFR', 'AMR', 'EAS', 'CSA', 'MID', 'TRANS']

# Load into a dictionary using a dict comprehension
pcs_dfs = {
    anc: pd.read_csv(f'pca/pca/{anc}_flashpca.pcs', sep='\t') 
    for anc in ancestries
}

# How to access a specific DataFrame later:
# df_within_pcs_EUR = pcs_dfs['EUR']

def process_covariate_and_save_data(df_case_control, df_within_pcs_ancestry, population):
    """
    Processes and saves covariate-related data for a specified ancestry group.

    Args:
        df_case_control (pd.DataFrame): DataFrame containing case-control data.
        df_within_pcs_ancestry (pd.DataFrame): DataFrame containing principal components
                                               and group information for a specific ancestry group.
        population (str): A string representing the ancestry population (e.g., 'EUR', 'AFR').
    """

    df_cohort_all_covariates = pd.merge(df_case_control, df_aou_v8_all_with_CNV[['IID', 'age', 'sex']], 
                                    how='left', left_on='IID', right_on='IID')
    
    # Create a list of IIDs from the ancestry PCs DataFrame
    ancestry_pc_iid_list = df_within_pcs_ancestry['IID'].tolist()
    
    # Filter the sliced DataFrame to include only individuals in the ancestry PC list
    df_ancestry = df_cohort_all_covariates[df_cohort_all_covariates['IID'].isin(ancestry_pc_iid_list)]
    
    # Merge the DataFrame with the PC data
    pc_columns = ['IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'group']
    df_ancestry = pd.merge(df_ancestry, df_within_pcs_ancestry[pc_columns], on='IID', how='left')
    
    # Keep only the "Unrelated_Inliers" group
    df_ancestry = df_ancestry[df_ancestry['group'] == 'Unrelated_Inliers']
    
    # Split data into cases (Case=1) and controls (Control=0)
    df_ancestry_controls = df_ancestry[df_ancestry['Control'] == 0]
    df_ancestry_cases    = df_ancestry[df_ancestry['Case'] == 1]
    
    # Drop the 'pheno' and 'group' columns from all DataFrames
    df_ancestry = df_ancestry.drop(columns=['pheno', 'group'])
    df_ancestry_controls = df_ancestry_controls.drop(columns=['pheno', 'group'])
    df_ancestry_cases = df_ancestry_cases.drop(columns=['pheno', 'group'])
    
    # Create the output directory if it doesn't exist
    output_dir = 'path/to/output/directory'  # Change this to your desired output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define file paths
    file_all      = os.path.join(output_dir, f'all_covariates_WITHIN_{population}_<phenotype>_age_sex_pcs.txt')
    file_controls = os.path.join(output_dir, f'covariates_WITHIN_{population}_<phenotype>_controls_age_sex_pcs.txt')
    file_cases    = os.path.join(output_dir, f'covariates_WITHIN_{population}_<phenotype>_cases_age_sex_pcs.txt')
    
    # Save the DataFrames to tab-separated files
    df_ancestry.to_csv(file_all, sep='\t', index=False)
    df_ancestry_controls.to_csv(file_controls, sep='\t', index=False)
    df_ancestry_cases.to_csv(file_cases, sep='\t', index=False)


populations = ['EUR', 'AFR', 'AMR', 'EAS', 'CSA', 'MID', 'TRANS']

# Assume df_case_control is already loaded
# df_case_control = pd.read_csv('case_control.txt', sep='\t')

for pop in populations:
    # Load the specific PC file for the current population
    # Adjust this line based on how your PC files are named
    # For example, if they are named like 'df_within_pcs_EUR.txt', 'df_within_pcs_AFR.txt', etc.
    df_within_pcs_ancestry = pd.read_csv(f'{pop}_ancestry.pcs', sep='\t')
    
    # Call the function for the current population
    process_covariate_and_save_data(df_case_control, df_within_pcs_ancestry, pop)

# check sample sizes for each covariate table in each population group
for pop in populations:
    df_all = pd.read_csv(f'path/to/output/directory/all_covariates_WITHIN_{pop}_<phenotype>_age_sex_pcs.txt', sep='\t')
    df_controls = pd.read_csv(f'path/to/output/directory/covariates_WITHIN_{pop}_<phenotype>_controls_age_sex_pcs.txt', sep='\t')
    df_cases = pd.read_csv(f'path/to/output/directory/covariates_WITHIN_{pop}_<phenotype>_cases_age_sex_pcs.txt', sep='\t')
    print(f"Population: {pop}")
    print(f"  All: {df_all.shape[0]}")
    print(f"  Controls: {df_controls.shape[0]}")
    print(f"  Cases: {df_cases.shape[0]}")

###################################################################################################################################
# create keep list of IIDs for samples for each ancestry group based on the covariate tables 
# and use that for downstream analyses (e.g., CNV burden analysis, CNV association analysis, etc.)
###################################################################################################################################
# code to create keep lists of IIDs for each ancestry group based on the covariate tables
for pop in populations:
    df_all = pd.read_csv(f'path/to/output/directory/all_covariates_WITHIN_{pop}_<phenotype>_age_sex_pcs.txt', sep='\t')
    keep_iids = df_all['IID'].tolist()
    with open(f'{pop}_keep_iids.txt', 'w') as f:
        for iid in keep_iids:
            f.write(f"{iid}\n")

#####################
# END
#####################