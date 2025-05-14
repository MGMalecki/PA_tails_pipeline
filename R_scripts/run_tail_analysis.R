# R_scripts/run_tail_analysis.R
# load libraries
library("stringr")
library("stringi")
library("tidyr")
library("dplyr")
source("R_scripts/run_AT_tail_analysis_df.R")
source("R_scripts/assign_type_and_class.R")
# load read1 and read2 samples
R1 <- read.table(snakemake@input[["tabR1"]], sep = "\t")
R2 <- read.table(snakemake@input[["tabR2"]], sep = "\t")
# strip UMI and safe it in new column
R2$UMI <- stringr::str_extract(R2$V1, "(?<=:)[^:]{15}$")
R2$V1 <- sub(":[^:]{15}$", "", R2$V1)
# merge read1 and read2 data
merged_dat <- dplyr::inner_join(
    dplyr::select(R1, V1, V2, V3, V4, V5, V7),
    dplyr::select(R2, V1, V2, V3, V4, V5, V6, UMI),
    by = "V1"
)
# assign column names
colnames(merged_dat)[1:11] <- c("read_name", "geneR1", "transR1", "flagR1","mapR1","CIGAR_R1",
                                "flagR2","chrR2","posR2","CIGAR_R2", "seq_R2")
# deduplicate the data
dedup_dat <- merged_dat %>% distinct(geneR1, UMI, .keep_all = TRUE)
dedup_dat$geneR1 <- sub("\\..*", "", dedup_dat$geneR1) # clean gene names
# calculate deduplication level
before_dedup <- merged_dat %>% count()
after_dedup <- dedup_dat %>% count()
dat <- as.data.frame(c("before_dedup" = before_dedup, "after_dedup" = after_dedup))
# merge data with reference
annotation <- readRDS("annotation/gene_annotations.rds")
dedup_dat <- dplyr::left_join(dedup_dat, annotation[,1:5], by = c("geneR1" = "gene_id"))
# analyse tail in the data
dedup_dat <- run_AT_tail_analysis_df(dedup_dat)
# assigna some qualifiers
dedup_dat <- assign_type_and_class(dedup_dat)
# load R2 reads that aligned end to end with bowtie
R2_2 <- read.table(snakemake@input[["tabR2_2"]], sep = "\t")
R2_2$V1 <- sub(":[^:]{15}$", "", R2_2$V1)
R2_2 <- R2_2 %>% filter(V2 == "unique")
grep_aligned <- dedup_dat %>% 
  filter(tail_len < 40) %>%
  filter(read_name %in% R2_2$V1 & type == "grep") %>%
  select(read_name)
# mark those reads in dedup dat
dedup_dat$al_ete_bt <- ifelse(dedup_dat$read_name %in% grep_aligned$read_name, "yes", "no")
# save output files
saveRDS(dedup_dat, snakemake@output[["big_table"]])
saveRDS(dedup_dat[,c(2,4,7:10,13:22,24:ncol(dedup_dat))], snakemake@output[["small_table"]])
write.csv(dat, snakemake@output[["deduplication"]])
