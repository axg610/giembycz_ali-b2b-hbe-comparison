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

```sh
module load R
R
```

```r
.libPaths("/home/alex.gao1/R")
setwd("/work/giembycz_lab/ALI-B2B-HBE_comparison/sleuth_ns/")

system("mkdir -p objects")

library(sleuth)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
```

```r prep sleuth inputs

t2g <- read_tsv("/work/giembycz_lab/ag_tools/Homo_sapiens.GRCh38.p13.cdna.all.070121.mart_export.txt") %>%
  select(target_id = 1, Gene = 2, gene_id = 4) %>% 
  distinct(target_id, gene_id, .keep_all = T) %>%
  na.omit()

meta = read_tsv("meta_ali-b2b-hbe-brushings_basal.txt") %>%
  filter(celltype != "brushings") %>%
  mutate(celltype = factor(celltype, levels = c("B2B", "HBE", "ALI"))) %>%

  # re-identify the technical sequencing runs since we scrubbed these out of the metafile for edgeR
  group_by(sample) %>%
  mutate(sample = paste0(sample, ".", row_number())) %>%
  ungroup()

new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

```

```r create sleuth object

so = sleuth_prep(
  sample_to_covariates = meta,
  target_mapping = t2g,
  full_model = ~celltype+time,
  gene_mode = TRUE,
  aggregation_column = "Gene",
  filter_fun = new_filter,
  num_cores = 4
)

saveRDS(so, "objects/so.rds")

```

```r fit models

so = sleuth_fit(
  so,
  ~celltype + time,
  method = "full"
  )

```


```r run tests with B2B as control

so <- sleuth_wt(
  so,
  which_beta = "celltypeHBE"
)

so <- sleuth_wt(
  so,
  which_beta = "celltypeALI"
)

res_HBE_vs_B2B <- sleuth_results(
  so,
  test = "celltypeHBE"
)

res_ALI_vs_B2B <- sleuth_results(
  so,
  test = "celltypeALI"
)

```

```r run tests with HBE as control

so$sample_to_covariates$celltype <- relevel(
  so$sample_to_covariates$celltype,
  ref = "HBE"
)

so <- sleuth_fit(
  so,
  ~celltype + time,
  method = "full"
)

so <- sleuth_wt(
  so,
  which_beta = "celltypeB2B"
)

so <- sleuth_wt(
  so,
  which_beta = "celltypeALI"
)

res_B2B_vs_HBE <- sleuth_results(
  so,
  test = "celltypeB2B"
)

res_ALI_vs_HBE <- sleuth_results(
  so,
  test = "celltypeALI"
)

```

```r run tests with ALI as control

so$sample_to_covariates$celltype <- relevel(
  so$sample_to_covariates$celltype,
  ref = "ALI"
)

so <- sleuth_fit(
  so,
  ~celltype + time,
  method = "full"
)

so <- sleuth_wt(
  so,
  which_beta = "celltypeB2B"
)

so <- sleuth_wt(
  so,
  which_beta = "celltypeHBE"
)

res_B2B_vs_ALI <- sleuth_results(
  so,
  test = "celltypeB2B"
)

res_HBE_vs_ALI <- sleuth_results(
  so,
  test = "celltypeHBE"
)

```



```r join results
all_celltype_results <- bind_rows(
  
  res_HBE_vs_B2B %>%
    mutate(contrast = "HBE_vs_B2B"),
  
  res_ALI_vs_B2B %>%
    mutate(contrast = "ALI_vs_B2B"),
  
  res_B2B_vs_HBE %>%
    mutate(contrast = "B2B_vs_HBE"),
  
  res_ALI_vs_HBE %>%
    mutate(contrast = "ALI_vs_HBE"),
  
  res_B2B_vs_ALI %>%
    mutate(contrast = "B2B_vs_ALI"),
  
  res_HBE_vs_ALI %>%
    mutate(contrast = "HBE_vs_ALI")
) %>%
  mutate(log2fold = b/log(2)) %>%
  select(Gene = target_id, contrast, log2fold, FDR = qval) %>%
  filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
  mutate(
    FDR = if_else(is.na(FDR), 1, FDR),
    log2fold = if_else(is.na(log2fold), 0, log2fold)
  ) %>%
  arrange(Gene, contrast) %>%
  distinct()

write_tsv(all_celltype_results, "basal_dea.txt")
```



