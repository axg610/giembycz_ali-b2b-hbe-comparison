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

files <- file.path(meta$path, "abundance.tsv")
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
# 3. run edgeR with TMM normalization, then pull cpms
# ------------------------------------------------------------

dge <- DGEList(counts = counts_bio)

dge <- calcNormFactors(dge, method = "TMM")

dge$samples

tmm_cpm <- cpm(dge, normalized.lib.sizes = TRUE)

# ------------------------------------------------------------
# 5. re-map gene ids and cleanup
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
    values_to = "tmm_count"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  relocate(tmm_count, .after = rep) %>%
  distinct()

saveRDS(tmm_df, "tmm_df.rds")

```












