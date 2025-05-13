# R_scripts/run_tail_analysis.R

library("stringr")
library("stringi")
library("tidyr")
library("dplyr")
source("R_scripts/run_AT_tail_analysis_df.R")
source("R_scripts/assign_type_and_class.R")

R1 <- read.table(snakemake@input[["tabR1"]], sep = "\t")
R2 <- read.table(snakemake@input[["tabR2"]], sep = "\t")

R2$UMI <- stringr::str_extract(R2$V1, "(?<=:)[^:]{15}$")
R2$V1 <- sub(":[^:]{15}$", "", R2$V1)

merged_dat <- dplyr::inner_join(
    dplyr::select(R1, V1, V2, V3, V4, V5, V7),
    dplyr::select(R2, V1, V2, V3, V4, V5, V6, UMI),
    by = "V1"
)

colnames(merged_dat)[1:11] <- c("read_name", "geneR1", "transR1", "flagR1","mapR1","CIGAR_R1",
                                "flagR2","chrR2","posR2","CIGAR_R2", "seq_R2")

dedup_dat <- merged_dat %>% distinct(geneR1, UMI, .keep_all = TRUE)
dedup_dat$geneR1 <- sub("\\..*", "", dedup_dat$geneR1)

annotation <- readRDS("annotation/gene_annotations.rds")
dedup_dat <- dplyr::left_join(dedup_dat, annotation[,1:5], by = c("geneR1" = "gene_id"))

dedup_dat <- run_AT_tail_analysis_df(dedup_dat)
dedup_dat <- assign_type_and_class(dedup_dat)

saveRDS(dedup_dat, snakemake@output[["big_table"]])
saveRDS(dedup_dat[,c(2,4,7:10,13:22,24,25)], snakemake@output[["small_table"]])

