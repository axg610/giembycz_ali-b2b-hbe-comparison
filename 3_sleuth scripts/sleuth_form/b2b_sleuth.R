.libPaths("/home/varuna.jayasinghe/R/x86_64-pc-linux-gnu-library/4.3")
setwd("/work/giembycz_lab/ALI-B2B-HBE_comparison/sleuth/")

library(dplyr)
library(readr)
library(sleuth)

#transcript ids to gene names
t2g <- read_tsv("/work/giembycz_lab/ag_tools/Homo_sapiens.GRCh38.p13.cdna.all.070121.mart_export.txt") %>%
  select(target_id = 1, Gene = 2, gene_id = 4) %>% 
  distinct(target_id, gene_id, .keep_all = T) %>%
  na.omit()

#sample file locations and associated covariates
s2c_2h <- read_tsv("b2b_meta.txt") %>% mutate(treatment = factor(treatment, levels = c("NS", "Form"))) %>% filter(time == "2h")
s2c_6h <- read_tsv("b2b_meta.txt") %>% mutate(treatment = factor(treatment, levels = c("NS", "Form"))) %>% filter(time == "6h")

#filter requiring at least 5 reads in at least 20% of samples
new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

#sleuth prep
so_2h <- sleuth_prep(sample_to_covariates = s2c_2h,
                     full_model = ~treatment,
                     target_mapping = t2g,
                     gene_mode = TRUE, 
                     aggregation_column = "Gene",
                     filter_fun = new_filter,
                     num_cores = 1)
so_2h <- sleuth_fit(so_2h, ~treatment)
so_2h <- sleuth_wt(so_2h, "treatmentForm")

so_6h <- sleuth_prep(sample_to_covariates = s2c_6h,
                     full_model = ~treatment,
                     target_mapping = t2g,
                     gene_mode = TRUE, 
                     aggregation_column = "Gene",
                     filter_fun = new_filter,
                     num_cores = 1)
so_6h <- sleuth_fit(so_6h, ~treatment)
so_6h <- sleuth_wt(so_6h, "treatmentForm")

#save sleuth objects
saveRDS(so_2h, "objects/b2b_2h.rds")
saveRDS(so_6h, "objects/b2b_6h.rds")

#assemble tpm table
tpm_2h <- kallisto_table(so_2h, use_filtered = FALSE) %>%
  select(Gene = target_id, sample, rep, time, treatment, tpm)

tpm_6h <- kallisto_table(so_6h, use_filtered = FALSE) %>%
  select(Gene = target_id, sample, rep, time, treatment, tpm)

tpm <- rbind(tpm_2h, tpm_6h)

write_tsv(tpm, "results/b2b_tpm.txt")

#assemble dea table
dea_2h <- sleuth_results(so_2h, "treatmentForm", "wt") %>%
  select(Gene = target_id, fold = b, FDR = qval) %>%
  mutate(fold = fold / log(2)) %>%
  mutate(celltype = "B2B", time = "2h", treatment = "Form", .before = 2)

dea_6h <- sleuth_results(so_6h, "treatmentForm", "wt") %>%
  select(Gene = target_id, fold = b, FDR = qval) %>%
  mutate(fold = fold / log(2)) %>%
  mutate(celltype = "B2B", time = "6h", treatment = "Form", .before = 2)

dea <- rbind(dea_2h, dea_6h)

write_tsv(dea, "results/b2b_dea.txt")