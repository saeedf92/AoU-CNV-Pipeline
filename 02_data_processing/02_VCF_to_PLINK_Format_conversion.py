import os
import pandas as pd
import subprocess
import time

# Set up the environment variables and directories
my_bucket = os.getenv('WORKSPACE_BUCKET')

SVvcfDIR     =  "gs://fc-aou-datasets-controlled/v8/wgs/short_read/structural_variants/vcf/full" # AoU CDRv8 SV VCFs
outputDIR    =  f"{my_bucket}/data/cnv_vcf_plink/plink_files"
logDIR       =  f"{my_bucket}/data/cnv_vcf_plink/logging"


# gathering the list of SV VCFs
command = "gsutil -u $GOOGLE_PROJECT ls " + SVvcfDIR + "/*.gz" 
vcf_list_gs = subprocess.check_output(command, shell=True)

decoded_vcf_list =  vcf_list_gs.decode('utf-8')
vcf_list_dir = decoded_vcf_list.strip().split('\n')

vcf_name = []
vcf_name_wo_ext = []

# looping through the list of VCF files and submitting a dsub job for each file to convert it to PLINK format
for vcf in vcf_list_dir:
    vcf_file_name = vcf.split('/')[-1]  # Extract the filename from the path
    vcf_name.append(vcf_file_name)      
    vcf_name_without_extension = vcf_file_name.replace(".vcf.gz", "")
    vcf_name_wo_ext.append(vcf_name_without_extension)
    
    # Building dsub command inside the loop to submit the jobs
    dsub_command = f'''
    source ~/aou_dsub.bash  

    aou_dsub \\
        --image biocontainer/plink2:alpha2.3_jan2020 \\
        --name "{vcf_file_name}" \\
        --boot-disk-size 300 \\
        --disk-size 300 \\
        --logging "{logDIR}" \\
        --input-recursive input="{SVvcfDIR}" \\
        --output-recursive output="{outputDIR}" \\
        --command 'plink2 \\
                    --vcf ${{input}}/{vcf_file_name} \\
                    --make-bed \\
                    --out ${{output}}/{vcf_name_without_extension}' 
    '''
    
    # Execute the command
    subprocess.run(dsub_command, shell=True, executable='/bin/bash')
    print(f'Job submitted for VCF file: {vcf_file_name}')
    time.sleep(1)