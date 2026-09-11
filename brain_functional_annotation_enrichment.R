#!/usr/bin/env Rscript

library(bigreadr)
library(tidyverse)

# Representative code for testing enrichment of significant WGS variants
# in brain ATAC-seq peaks.
#
# Project-specific input paths have been replaced with placeholders.
# This script is not intended to serve as a standalone executable pipeline.


# ==============================================================================
# Input data
# ==============================================================================

# Variant IDs included in the REGENIE variant-level testing universe.
all_variants <- fread2('PATH_TO_ALL_TESTED_REGENIE_VARIANTS')$ID

# Brain ATAC-seq annotation file containing:
#   V1: variant ID
#   V3: ATAC-seq peak category ("neuron", "glia", or "all")
atac <- fread2("PATH_TO_BRAIN_ATAC_SEQ_PEAK_ANNOTATIONS")

neuron_variants <- unique(atac$V1[atac$V3 == 'neuron'])         # ATAC seq peaks in neuron
glia_variants <- unique(atac$V1[atac$V3 == 'glia'])             # ATAC seq peaks in glia
all_cell_types_variants <- unique(atac$V1[atac$V3 == 'all'])    # ATAC seq peaks in all cell types


# ==============================================================================
# Enrichment test
# ==============================================================================

fisher_enrichment_test <- function(significant_variants, nonsignificant_variants, annotation_variants) {
      
  a <- sum(significant_variants %in% annotation_variants)                                   
  b <- length(significant_variants) - a      
  c <- sum(nonsignificant_variants %in% annotation_variants)                                
  d <- length(nonsignificant_variants) - c
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  colnames(tab) <- c("in_group", "not_in_group")
  rownames(tab) <- c("significant", "not_significant")
  f <- fisher.test(tab, alternative = "greater")

  observed_prop <- a / (a + b)
  background_prop <- (a + c) / (a + b + c + d)

  return(list(
    sig_in = a,
    sig_not_in = b,
    not_sig_in = c,
    not_sig_not_in = d,
    odds_ratio = unname(f$estimate),
    p_value = f$p.value,
    conf_low = f$conf.int[1],
    conf_high = f$conf.int[2],
    observed_prop = observed_prop,
    background_prop = background_prop,
    fold_enrichment = observed_prop / background_prop
  ))
}


# ==============================================================================
# Enrichment analyses
# ==============================================================================

# significant variant-level associations for regional brain volumes, DTI parameters, and ICA nodes
result_files <- c(
  "Regional brain volumes" = "PATH_TO_SIGNIFICANT_VARIANT_LEVEL_REGIONAL_BRAIN_VOLUMES_RESULTS",
  "DTI parameters" = "PATH_TO_SIGNIFICANT_VARIANT_LEVEL_DTI_PARAMETERS_RESULTS",
  "ICA nodes" = "PATH_TO_SIGNIFICANT_VARIANT_LEVEL_ICA_NODES_RESULTS"
)

annotation_sets <- list(
  neuron = neuron_variants,
  glia = glia_variants,
  all_cell_types = all_cell_types_variants
)

enrichment_results <- data.frame()

for (imaging_modality in names(result_files)) {
  association_results <- fread2(result_files[[imaging_modality]])
  
  for (trait in unique(association_results$trait)) {
    significant_variants <- association_results$ID[association_results$trait == trait]
    nonsignificant_variants <- all_variants[!all_variants %in% significant_variants]
    
    for (cell_type in names(annotation_sets)) {
      result <- fisher_enrichment_test(significant_variants, nonsignificant_variants, annotation_sets[[cell_type]])
      result$trait <- trait
      result$cell_type <- cell_type
      result$imaging_modality <- imaging_modality
      
      enrichment_results <- bind_rows(enrichment_results, result)
    }
  }
}

# Multiple testing
enrichment_results$p_bonferroni <- p.adjust(enrichment_results$p_value, method='bonferroni')
enrichment_results$significant <- enrichment_results$p_bonferroni < 0.05

write.csv(enrichment_results, 'brain_ATAC_seq_peak_enrichment_results.csv', row.names = FALSE)
