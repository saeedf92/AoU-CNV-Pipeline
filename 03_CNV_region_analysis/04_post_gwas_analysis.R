library(tidyverse)
library(qqman)

# Set working directory to the location of the results
setwd("/path/to/your/analysis/directory")

# concatinate gwas results
read_and_combine_firth <- function(input_path, input_pattern = "*.glm.firth") {
    
    file_list <- list.files(path = input_path, pattern = input_pattern, full.names = TRUE)

    read_firth <- function(file) { 
      read.csv(file, header = TRUE, sep = "\t")
    }

    combined_data <- do.call(rbind, lapply(file_list, read_firth))
    return(combined_data)
}

ancestries <- c("EUR", "AFR", "AMR", "EAS", "CAS", "TRANS" )
combined_data_list <- list()
for (anc in ancestries) {
  path <- paste0("./AUD_within_", anc, "_Results/")
  combined_data_list[[paste0("combined_data_", anc)]] <- read_and_combine_firth(
    input_path = path, 
    input_pattern = "*glm.firth"
  )
}
# How to access your data now:
# combined_data_list$combined_data_EUR
# combined_data_list$combined_data_AFR
# combined_data_list$combined_data_AMR
# combined_data_list$combined_data_EAS
# combined_data_list$combined_data_CAS
# combined_data_list$combined_data_TRANS


firth_DelDup_modif <- combined_data_<ANCESTRY> %>%
    filter(`TEST` == "ADD") %>%
    filter(`ALT` == '<DUP>' | `ALT` == '<DEL>' ) %>%
    filter(`X.CHROM` != "X") %>%
    filter(!is.na(`P`))

firth_DelDup_modif$X.CHROM <- as.integer(firth_DelDup_modif$X.CHROM)

firth_DelDup_modif <- firth_DelDup_modif %>%
    rename(CHR = X.CHROM)

# save sum stat
write.table(firth_DelDup_modif, 
          file = "./<PHENOTYPE>_<ANCESTRY>_SumStat.tsv",
          sep = "\t",
          row.names = FALSE,
          quote = FALSE)

is.numeric(firth_DelDup_modif$CHR)

# QQ plot
qq(firth_DelDup_modif$P)

# Manhattan plot
# Create a temporary numeric version of the 'P' column
temp_p <- as.numeric(as.character(firth_DelDup_modif$P))
# Find the rows where the conversion failed (i.e., became NA)
problem_rows <- firth_DelDup_modif[is.na(temp_p), ]
# Look at the values in the 'P' column from these rows
# This will show you the exact strings causing the issue
# print(problem_rows$P)

png("./manhattan_plot_<PHENOTYPE>_<ANCESTRY>_firth.png", width=1000, height=600) 
manhattan(firth_DelDup_modif, chr="CHR", bp="POS", snp="ID", p="P", main="CNVs Firth Regression", 
          annotatePval = 0.00001, col = c("#003049", "#c1121f"),
          cex.axis = 0.8, las=0)
dev.off()

# FDR adjustment
firth_DelDup_modif <- firth_DelDup_modif %>%
    mutate(FDR = p.adjust((P), method = "BH")) %>%
    mutate(BONF = p.adjust((P), method = "bonferroni")) %>%
    arrange(FDR)
# head(firth_DelDup_modif, 5)

# Lambda and lambda1000
z_stats <- firth_DelDup_modif$Z_STAT
chi_squared <- z_stats^2
lambda     <-  median(chi_squared)/qchisq(0.5,1)
lambda1000 <-  1 + (lambda - 1) * (1/n_cases + 1/n_controls) / (1/1000 + 1/1000)

# Annotating CNVs
cnv_annotation_ref <- read.csv("annotation/AoU_srWGS_SV.v8.sites_only.vcf_Extracted_ID_POS_SVLEN_AF_cytoband_genes.tsv", 
                               header = TRUE,  sep = "\t")
annotated_df <- firth_DelDup_modif %>%
  left_join(cnv_annotation_ref %>% select(ID, start, end, `length..bp.`,  MAF, cytoband, gene), by = "ID")
