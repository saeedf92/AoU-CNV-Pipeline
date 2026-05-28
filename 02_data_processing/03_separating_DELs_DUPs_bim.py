import os
import pandas as pd

bucket = os.getenv('WORKSPACE_BUCKET')

# look one of .bim files head
# ! gsutil cat {bucket}/data/cnv_vcf_plink/plink_files/AoU_srWGS_SV.v8.chr1.bim | head -n 20

# copy .bim files to local
# ! gsutil -m cp {bucket}/data/cnv_vcf_plink/plink_files/*.bim ./aud/bim_files

# read .bim files, separate DELs and DUPs, and save them as separate files
bim_directory = "aud/bim_files/"

for filename in os.listdir(bim_directory):
    if filename.endswith(".bim"):
        file_path = os.path.join(bim_directory, filename)

        df = pd.read_csv(file_path, delim_whitespace=True, header=None)

        # Filter rows where the 5th column contains <DEL> or <DUP>
        dels = df[df[4] == "<DEL>"]
        dups = df[df[4] == "<DUP>"]
        
        base_filename = os.path.splitext(filename)[0]
        
        dels_file_path = os.path.join(bim_directory, f"{filename}_DELs.txt")
        dups_file_path = os.path.join(bim_directory, f"{filename}_DUPs.txt")

        dels.to_csv(dels_file_path, sep='\t', header=False, index=False)
        dups.to_csv(dups_file_path, sep='\t', header=False, index=False)

        print(f"Processed {filename}: DELs saved to {dels_file_path}, DUPs saved to {dups_file_path}")


# After running the above code, you will have separate files for DELs and DUPs in the same directory. 
# You can then copy these files back to the bucket for downstream analysis.
# ! gsutil cp aud/bim_files/*.txt {bucket}/plink_files/