Hi! Alex.

This is what we require (I think):

1. Comparative resting (unstimulated) expression of genes in B2B, HBE, ALIs and bronchial brushings expressed as TMMs.
 

2. DEA of resting gene expression in HBEs vs. ALIs, and HBEs vs. B2Bs (I want to generate volcano plots etc to illustrate, at a gene level, that these cells have basal cell vs secretory/ciliated cell characteristics).
 
Mark

# setup arc

```sh
ssh alex.gao1@arc.ucalgary.ca
salloc -c 8 --mem 128GB --time 05:00:00
```

# compare resting expression via TMM adjustment

```sh
module load R
R
```

```r
.libPaths("/home/alex.gao1/R")
setwd("/work/giembycz_lab/ALI-B2B-HBE_comparison/basal_tmm_adjustment/")

library(tximport)
library(edgeR)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)

meta = read_tsv("meta_ali-b2b-hbe-brushings_basal.txt")

tx2gene <- read_tsv("/work/newton_lab/mm_analysis/ref_sequence/GRCh38.p13_hgnc_t2g.txt") %>%
    select(TXNAME = target_id, GENEID = gene_id)

files <- file.path(meta$path, "abundance.h5")
names(files) <- meta$sample

# ------------------------------------------------------------
# 1. Run tximport
# ------------------------------------------------------------

txi <- tximport(
  files = files,
  type = "kallisto",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = TRUE
)

# ------------------------------------------------------------
# 2. collapse technical sequencing runs to biological samples
# ------------------------------------------------------------

# Biological sample IDs
sample_id <- meta$sample

# Collapse counts by biological sample
counts_bio <- sapply(unique(sample_id), function(s) {
  idx <- sample_id == s
  rowSums(txi$counts[, idx, drop = FALSE])
})

# Preserve gene IDs
rownames(counts_bio) <- rownames(txi$counts)








# ------------------------------------------------------------
# 3b. TMM: run edgeR with TMM normalization, then pull cpms
# ------------------------------------------------------------

dge <- DGEList(counts = counts_bio)

dge <- calcNormFactors(dge, method = "TMM")

write_tsv(
  dge$samples %>% rownames_to_column("sample"), 
  "tmm_normfactors.txt"
  )

tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

# ------------------------------------------------------------
# 4b. re-map gene ids and cleanup
# ------------------------------------------------------------

gene_map <- read_tsv(
  "/work/newton_lab/mm_analysis/ref_sequence/GRCh38.p13_hgnc_t2g.txt"
) %>%
  select(GENEID = gene_id, Gene = gene_symbol) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(GENEID, Gene) %>%
  group_by(GENEID) %>%
  filter(n() == 1) %>%
  ungroup()

sample_info <- meta %>%
  distinct(sample, celltype, treatment, time, rep)

tmm_df <- as.data.frame(tmm_cpm) %>%
  rownames_to_column("GENEID") %>%
  left_join(gene_map, by = "GENEID") %>%
  filter(!is.na(Gene)) %>%
  relocate(Gene, .before = 1) %>%
  select(-GENEID) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  relocate(value, .after = rep) %>%
  distinct()

saveRDS(tmm_df, "tmm_df.rds")




# ------------------------------------------------------------
# 3c. GeTMM: normalize to effective length, then run TMM
# ------------------------------------------------------------

# gene lengths corresponding to each biological sample
length_bio <- sapply(unique(sample_id), function(s) {
    idx <- sample_id == s
    rowMeans(txi$length[, idx, drop = FALSE])
})

# preserve gene id's
rownames(length_bio) <- rownames(txi$length)

# correct counts_bio with gene lengths
counts_length_corrected <- counts_bio / length_bio


# run edgeR with TMM normalization

dge_length_corrected = DGEList(counts = counts_length_corrected)

dge_length_corrected <- calcNormFactors(
    dge_length_corrected,
    method = "TMM"
)

write_tsv(
  dge_length_corrected$samples %>% rownames_to_column("sample"), 
  "geTMM_normfactors.txt"
  )

geTMM_cpm = cpm(
  dge_length_corrected,
  normalized.lib.sizes = TRUE
)

# ------------------------------------------------------------
# 4c. re-map gene ids and cleanup
# ------------------------------------------------------------

gene_map <- read_tsv(
  "/work/newton_lab/mm_analysis/ref_sequence/GRCh38.p13_hgnc_t2g.txt"
) %>%
  select(GENEID = gene_id, Gene = gene_symbol) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(GENEID, Gene) %>%
  group_by(GENEID) %>%
  filter(n() == 1) %>%
  ungroup()

sample_info <- meta %>%
  distinct(sample, celltype, treatment, time, rep)

geTMM_df <- as.data.frame(geTMM_cpm) %>%
  rownames_to_column("GENEID") %>%
  left_join(gene_map, by = "GENEID") %>%
  filter(!is.na(Gene)) %>%
  relocate(Gene, .before = 1) %>%
  select(-GENEID) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  relocate(value, .after = rep) %>%
  distinct()

saveRDS(geTMM_df, "geTMM_df.rds")



# ------------------------------------------------------------
# 3d. TPM
# ------------------------------------------------------------

tpm_bio <- sapply(unique(sample_id), function(s) {
  idx <- sample_id == s
  
  # TPM/abundance is not additive across technical runs,
  # so take the mean across runs
  rowMeans(txi$abundance[, idx, drop = FALSE])
})

rownames(tpm_bio) <- rownames(txi$abundance)

tpm_df <- as.data.frame(tpm_bio) %>%
  rownames_to_column("GENEID") %>%
  left_join(gene_map, by = "GENEID") %>%
  filter(!is.na(Gene)) %>%
  relocate(Gene, .before = 1) %>%
  select(-GENEID) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "sample",
    values_to = "value"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  relocate(value, .after = rep) %>%
  distinct()

saveRDS(tpm_df, "tpm_df.rds")
```












