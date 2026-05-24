"""
AoU Universal Phenotype Cohort Extraction Framework
Author: [Saeed Farajzadeh Valilou]
Description: A modular, parameterized template to programmatically extract 
             any target cohort from the All of Us Controlled/Registered Tier 
             using condition codes, drug exposures, and survey parameters.
             It is 
"""

import os
import sys
import logging
import pandas as pd
from google.cloud import bigquery

# Configure terminal logging (optional)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ==============================================================================
# USER CONFIGURATION BLOCK: DEFINE YOUR PHENOTYPE TARGETS HERE
# ==============================================================================
# Example settings below are set to placeholders. Update these for your trait.

# 1. Condition Codes: Add ICD-9 and ICD-10 prefixes (e.g., 'F10' for AUD) or exact OMOP concept IDs
CONDITION_ICD_PREFIXES = ['%REPLACE_WITH_ICD_PREFIX%']  # e.g., '%F10%', '%E11%'
CONDITION_CONCEPT_IDS  = [000000, 111111]               # e.g., Standard OMOP IDs

# 2. Medications: Add drug names in uppercase (e.g., 'Neltrexone') or use concept IDs if preferred
MEDICATION_KEYWORDS   = ['%REPLACE_WITH_DRUG_NAME%', '%ANOTHER_DRUG%']

# 3. Survey Data: Add uppercase keywords found in questionnaire maps 
SURVEY_KEYWORDS       = ['%REPLACE_WITH_SURVEY_KEYWORD%', '%LIFESTYLE%']
# ==============================================================================

def build_universal_query(cdr_dataset_id: str) -> str:
    """Dynamically generates the SQL framework injection based on user parameters."""
    
    # Format list elements securely for SQL parsing
    icd_clauses = " OR ".join([f"c.concept_code LIKE '{x}'" for x in CONDITION_ICD_PREFIXES])
    concept_id_str = ", ".join([str(x) for x in CONDITION_CONCEPT_IDS]) or "0"
    
    med_clauses = " OR ".join([f"UPPER(c.concept_name) LIKE '{x}'" for x in MEDICATION_KEYWORDS])
    survey_clauses = " OR ".join([f"UPPER(c.concept_name) LIKE '{x}'" for x in SURVEY_KEYWORDS])

    return f"""
    WITH target_conditions AS (
        SELECT DISTINCT 
            co.person_id, 
            'Condition' AS data_source, 
            co.condition_start_date AS event_date,
            COALESCE(c.concept_code, CAST(co.condition_concept_id AS STRING)) AS diagnostic_code
        FROM `{cdr_dataset_id}.condition_occurrence` co
        JOIN `{cdr_dataset_id}.concept` c ON co.condition_concept_id = c.concept_id
        WHERE ({icd_clauses})
           OR co.condition_concept_id IN ({concept_id_str}) 
    ),
    
    target_medications AS (
        SELECT DISTINCT 
            de.person_id, 
            'Medication' AS data_source, 
            de.drug_exposure_start_date AS event_date,
            COALESCE(c.concept_code, CAST(de.drug_concept_id AS STRING)) AS diagnostic_code
        FROM `{cdr_dataset_id}.drug_exposure` de
        JOIN `{cdr_dataset_id}.concept` c ON de.drug_concept_id = c.concept_id
        WHERE ({med_clauses})
    ),
    
    target_surveys AS (
        SELECT DISTINCT 
            o.person_id, 
            'Survey' AS data_source, 
            o.observation_date AS event_date,
            COALESCE(c.concept_code, CAST(o.observation_concept_id AS STRING)) AS diagnostic_code
        FROM `{cdr_dataset_id}.observation` o
        JOIN `{cdr_dataset_id}.concept` c ON o.observation_concept_id = c.concept_id
        WHERE c.concept_class_id = 'Question'
          AND ({survey_clauses})
    ),
    
    combined_cohort AS (
        SELECT person_id, data_source, event_date, diagnostic_code FROM target_conditions
        UNION DISTINCT
        SELECT person_id, data_source, event_date, diagnostic_code FROM target_medications
        UNION DISTINCT
        SELECT person_id, data_source, event_date, diagnostic_code FROM target_surveys
    ),
    
    aggregated_events AS (
        SELECT 
            person_id, 
            MIN(event_date) AS earliest_phenotype_date,
            STRING_AGG(DISTINCT diagnostic_code, ', ' ORDER BY diagnostic_code) AS matched_codes
        FROM combined_cohort
        GROUP BY person_id
    )
    
    SELECT 
        ae.person_id,
        ae.earliest_phenotype_date,
        ae.matched_codes,
        
        -- Standardized Core Demographics
        p.gender_concept_id,
        COALESCE(gc.concept_name, 'Unknown/Not Reported') AS gender,
        
        EXTRACT(DATE FROM p.birth_datetime) AS date_of_birth,
        
        p.race_concept_id,
        COALESCE(rc.concept_name, 'Unknown/Not Reported') AS race,
        
        p.ethnicity_concept_id,
        COALESCE(ec.concept_name, 'Unknown/Not Reported') AS ethnicity,
        
        p.sex_at_birth_concept_id,
        COALESCE(sc.concept_name, 'Unknown/Not Reported') AS sex_at_birth
        
    FROM aggregated_events ae
    JOIN `{cdr_dataset_id}.person` p ON ae.person_id = p.person_id
    LEFT JOIN `{cdr_dataset_id}.concept` gc ON p.gender_concept_id = gc.concept_id
    LEFT JOIN `{cdr_dataset_id}.concept` rc ON p.race_concept_id = rc.concept_id
    LEFT JOIN `{cdr_dataset_id}.concept` ec ON p.ethnicity_concept_id = ec.concept_id
    LEFT JOIN `{cdr_dataset_id}.concept` sc ON p.sex_at_birth_concept_id = sc.concept_id
    """

def extract_phenotype_cohort(output_path: str):
    """Executes the custom query generation and downloads data to a local TSV."""
    try:
        billing_project_id = os.environ['GOOGLE_PROJECT']
        cdr_dataset_id = os.environ['WORKSPACE_CDR']
    except KeyError as e:
        logging.error(f"Missing environment variable: {e}. Run inside an active AoU Workbench Session.")
        sys.exit(1)
        
    client = bigquery.Client(project=billing_project_id)
    query = build_universal_query(cdr_dataset_id)
    
    logging.info("Running universal multi-domain cohort build query...")
    try:
        cohort_df = client.query(query).to_dataframe()
    except Exception as e:
        logging.error(f"BigQuery generation error: {e}")
        sys.exit(1)
        
    logging.info(f"Complete! Extracted {len(cohort_df)} unique case profiles.")
    cohort_df.to_csv(output_path, sep="\t", index=False)
    logging.info(f"Matrix file successfully outputted to: '{output_path}'")

if __name__ == "__main__":
    target_output = sys.argv[1] if len(sys.argv) > 1 else "custom_phenotype_cohort.tsv"
    extract_phenotype_cohort(output_path=target_output)