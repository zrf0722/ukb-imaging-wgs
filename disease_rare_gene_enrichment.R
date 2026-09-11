#!/usr/bin/env Rscript

library(bigreadr)
library(tidyverse)

# Representative code for testing enrichment of significant and conditionally
# significant gene-based associations in disease-associated rare-variant genes.
#
# Project-specific input paths have been replaced with placeholders.
# This script is not intended to serve as a standalone executable pipeline.


# Significant gene-based associations for regional brain volumes, DTI parameters, and ICA nodes
regional_brain_volumes <- fread2("PATH_TO_SIGNIFICANT_GENE_BASED_REGIONAL_BRAIN_VOLUMES_RESULTS")
dti_parameters <- fread2("PATH_TO_SIGNIFICANT_GENE_BASED_DTI_PARAMETERS_RESULTS")
ica_nodes <- fread2("PATH_TO_SIGNIFICANT_GENE_BASED_ICA_NODES_RESULTS")

# Gene-based associations that remained significant after conditioning on independent significant 
# GWAS variants and had no significant GWAS variant within 100 kb of the gene
regional_brain_volumes_cond <- fread2("PATH_TO_CONDITIONALLY_SIGNIFICANT_GENE_BASED_REGIONAL_BRAIN_VOLUMES_RESULTS")
dti_parameters_cond <- fread2("PATH_TO_CONDITIONALLY_SIGNIFICANT_GENE_BASED_DTI_PARAMETERS_RESULTS")
ica_nodes_cond <- fread2("PATH_TO_CONDITIONALLY_SIGNIFICANT_GENE_BASED_ICA_NODES_RESULTS")

disease_gene_enrichment <- function(modality_res, disease) {
  all_traits <- unique(modality_res$trait)
  res_all <- data.frame()
  
  if (disease == 'autism_spectrum_disorder') {
    # 436 high-confidence autism spectrum disorder genes by SFARI: https://journals.biologists.com/dmm/article/3/3-4/133/2349/SFARI-Gene-an-evolving-database-for-the-autism
    asd <- fread2('SFARI-Gene_genes_05-01-2026release_05-16-2026export.csv')
    sfari_1_or_syndromic <- asd %>%
      filter(`gene-score` == 1 | syndromic == 1) %>%
      distinct(`gene-symbol`, .keep_all = TRUE)
    disease_genes <- sfari_1_or_syndromic$`gene-symbol`
  } else if (disease == 'developmental_disorder') {
    # 285 developmental disorder associated genes by Kaplanis et al: https://www.nature.com/articles/s41586-020-2832-5
    developmental_disorder_genes <- c('ABL1', 'ADCY5', 'ADNP', 'AFF3', 'AGO1', 'AHDC1', 'ALG13', 'AMER1', 'ANKRD11', 'AP2S1', 'ARF1', 'ARHGAP35', 'ARID1A', 'ARID1B', 'ARID2', 'ASXL1', 'ASXL3', 'ATP1A3', 'ATP6V0A1', 'ATP6V1A', 'ATRX', 'AUTS2', 'BCL11A', 'BCL11B', 'BCOR', 'BRAF', 'BRPF1', 'CACNA1A', 'CACNA1C', 'CACNA1E', 'CACNA1G', 'CAMK2B', 'CAMK2G', 'CAMTA1', 'CASK', 'CBL', 'CDK13', 'CDK8', 'CDKL5', 'CHAMP1', 'CHD2', 'CHD3', 'CHD4', 'CHD7', 'CHD8', 'CLCN4', 'CLTC', 'CNKSR2', 'CNOT3', 'COG4', 'COL2A1', 'COL4A1', 'COL4A3BP', 'COL6A1', 'CREBBP', 'CSNK2A1', 'CSNK2B', 'CTBP1', 'CTCF', 'CTNNB1', 'DDX23', 'DDX3X', 'DEAF1', 'DHDDS', 'DLG4', 'DNM1', 'DNM1L', 'DNMT3A', 'DPF2', 'DSP', 'DYNC1H1', 'DYRK1A', 'EBF3', 'EEF1A2', 'EFTUD2', 'EHMT1', 'EP300', 'ERF', 'FBN1', 'FBXO11', 'FBXW7', 'FGF12', 'FGFR2', 'FGFR3', 'FOXG1', 'FOXP1', 'FOXP2', 'GABRA1', 'GABRB2', 'GABRB3', 'GATA3', 'GATAD2B', 'GIGYF1', 'GNAI1', 'GNAO1', 'GNAS', 'GNB1', 'GNB2', 'GRIA2', 'GRIN1', 'GRIN2A', 'GRIN2B', 'H3F3A', 'HDAC4', 'HDAC8', 'HECW2', 'HIST1H1E', 'HIVEP2', 'HK1', 'HNRNPD', 'HNRNPK', 'HNRNPU', 'HRAS', 'HUWE1', 'IQSEC2', 'IRF2BPL', 'ITPR1', 'KANSL1', 'KAT6A', 'KAT6B', 'KCNA2', 'KCNB1', 'KCND3', 'KCNH1', 'KCNK3', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KDM5B', 'KDM5C', 'KDM6A', 'KDM6B', 'KIDINS220', 'KIF11', 'KIF1A', 'KIF5C', 'KLF7', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KMT2E', 'KMT5B', 'KRAS', 'LZTR1', 'MAGEL2', 'MAP2K1', 'MAP3K7', 'MAST1', 'MBD5', 'MECP2', 'MED12', 'MED13L', 'MEF2C', 'MEIS2', 'MFN2', 'MIB1', 'MMGT1', 'MN1', 'MORC2', 'MSL2', 'MSL3', 'MTOR', 'MYCN', 'MYT1L', 'NAA10', 'NAA15', 'NACC1', 'NEXMIF', 'NFIA', 'NFIX', 'NOTCH1', 'NR2F1', 'NR4A2', 'NSD1', 'NSD2', 'ODC1', 'OGT', 'PACS1', 'PACS2', 'PBX1', 'PCDH19', 'PDHA1', 'PHF21A', 'PHF6', 'PHIP', 'PIK3CA', 'PIK3R1', 'PIK3R2', 'POGZ', 'POU3F3', 'PPM1D', 'PPP1CB', 'PPP2CA', 'PPP2R1A', 'PPP2R5D', 'PRKAR1B', 'PRPF8', 'PRR12', 'PSMC5', 'PTCH1', 'PTEN', 'PTPN11', 'PUF60', 'PURA', 'QRICH1', 'RAB14', 'RAC1', 'RAI1', 'RERE', 'RHOBTB2', 'RIT1', 'RPS6KA3', 'SAMD9', 'SATB1', 'SATB2', 'SCN1A', 'SCN2A', 'SCN8A', 'SET', 'SETBP1', 'SETD1A', 'SETD2', 'SETD5', 'SHANK3', 'SHOC2', 'SIN3A', 'SKI', 'SLC2A1', 'SLC35A2', 'SLC6A1', 'SMAD4', 'SMARCA2', 'SMARCA4', 'SMARCB1', 'SMC1A', 'SNAP25', 'SON', 'SOX10', 'SOX11', 'SOX5', 'SPAST', 'SPEN', 'SPTAN1', 'SRCAP', 'SRRM2', 'STAG2', 'STXBP1', 'SYNGAP1', 'SYT1', 'TAB2', 'TAF1', 'TAOK1', 'TBL1XR1', 'TCF20', 'TCF4', 'TCF7L2', 'TFE3', 'TLK2', 'TMEM106B', 'TRAF7', 'TRIO', 'TRIP12', 'TRPM3', 'TRRAP', 'U2AF2', 'UBE3A', 'UBTF', 'UPF1', 'USP7', 'USP9X', 'VAMP2', 'WAC', 'WDR26', 'WDR37', 'WDR45', 'YY1', 'ZBTB18', 'ZBTB20', 'ZC4H2', 'ZEB2', 'ZFHX4', 'ZMYND11', 'ZNF148', 'ZNF292')
    disease_genes <- developmental_disorder_genes
  } 

  for (trait in all_traits) {
    res_tmp <- modality_res[modality_res$trait == trait, ]
    genes_tmp <- unique(res_tmp$gene)

    ## One-sided hypergeometric enrichment test
    N <- 18371                                  # genes included in the REGENIE gene-based testing universe
    K <- length(genes_tmp)                      # significant genes
    n <- length(disease_genes)                  # disease genes
    x <- sum(genes_tmp %in% disease_genes)      # overlap

    a <- x              # disease and significant
    b <- n - x          # disease and not significant
    c <- K - x          # non-disease and significant
    d <- N - n - c      # non-disease and not significant
  
    odds_ratio <- (a * d) / (b * c)  # Odds ratio
    expected <- n * K / N # Expected overlap
    fold_enrichment <- x / expected # Fold enrichment
    p_hyper <- phyper(x - 1, K, N - K, n, lower.tail = FALSE) 
    overlap_genes <- genes_tmp[genes_tmp %in% disease_genes]
    
    res_all <- rbind.data.frame(
      res_all,
      data.frame(
        trait = trait,
        n_genes = length(genes_tmp),
        n_overlap = x,
        a = a,
        b = b,
        c = c,
        d = d,
        odds_ratio = odds_ratio,
        expected = expected,
        fold_enrichment = fold_enrichment,
        p_hyper = p_hyper,
        overlapped_genes = paste(overlap_genes, collapse = ","),
        stringsAsFactors = FALSE
      )
    )
  }
  return(res_all)
}

##### Autism spectrum disorder 
# All significant
regional_brain_volumes_enrichment_asd <- disease_gene_enrichment(regional_brain_volumes, 'autism_spectrum_disorder')
dti_parameters_enrichment_asd <- disease_gene_enrichment(dti_parameters, 'autism_spectrum_disorder')
ica_nodes_enrichment_asd <- disease_gene_enrichment(ica_nodes, 'autism_spectrum_disorder')

# Conditionally significant and do not have genome-wide significant GWAS variants within 100kb window
regional_brain_volumes_enrichment_cond_asd <- disease_gene_enrichment(regional_brain_volumes_cond, 'autism_spectrum_disorder')
dti_parameters_enrichment_cond_asd <- disease_gene_enrichment(dti_parameters_cond, 'autism_spectrum_disorder')
ica_nodes_enrichment_cond_asd <- disease_gene_enrichment(ica_nodes_cond, 'autism_spectrum_disorder')

##### Developmental disorder
# All significant
regional_brain_volumes_enrichment_dd <- disease_gene_enrichment(regional_brain_volumes, 'developmental_disorder')
dti_parameters_enrichment_dd <- disease_gene_enrichment(dti_parameters, 'developmental_disorder')
ica_nodes_enrichment_dd <- disease_gene_enrichment(ica_nodes, 'developmental_disorder')

# Conditionally significant and do not have genome-wide significant GWAS variants within 100kb window
regional_brain_volumes_enrichment_cond_dd <- disease_gene_enrichment(regional_brain_volumes_cond, 'developmental_disorder')
dti_parameters_enrichment_cond_dd <- disease_gene_enrichment(dti_parameters_cond, 'developmental_disorder')
ica_nodes_enrichment_cond_dd <- disease_gene_enrichment(ica_nodes_cond, 'developmental_disorder')

# Add imaging modality
regional_brain_volumes_enrichment_asd$imaging_modality <- regional_brain_volumes_enrichment_cond_asd$imaging_modality <- 
  regional_brain_volumes_enrichment_dd$imaging_modality <- regional_brain_volumes_enrichment_cond_dd$imaging_modality <- 'Regional brain volumes'
dti_parameters_enrichment_asd$imaging_modality <- dti_parameters_enrichment_cond_asd$imaging_modality <- 
  dti_parameters_enrichment_dd$imaging_modality <- dti_parameters_enrichment_cond_dd$imaging_modality <- 'DTI parameters'
ica_nodes_enrichment_asd$imaging_modality <- ica_nodes_enrichment_cond_asd$imaging_modality <- 
  ica_nodes_enrichment_dd$imaging_modality <- ica_nodes_enrichment_cond_dd$imaging_modality <- 'ICA nodes'

# Add disease
regional_brain_volumes_enrichment_asd$disease <- regional_brain_volumes_enrichment_cond_asd$disease <- dti_parameters_enrichment_asd$disease <- 
  dti_parameters_enrichment_cond_asd$disease <- ica_nodes_enrichment_asd$disease <- ica_nodes_enrichment_cond_asd$disease <- 'Autism spectrum disorder'
regional_brain_volumes_enrichment_dd$disease <- regional_brain_volumes_enrichment_cond_dd$disease <- dti_parameters_enrichment_dd$disease <- 
  dti_parameters_enrichment_cond_dd$disease <- ica_nodes_enrichment_dd$disease <- ica_nodes_enrichment_cond_dd$disease <- 'Developmental disorder'

# Add significant or conditionally significant info
regional_brain_volumes_enrichment_asd$group <- regional_brain_volumes_enrichment_dd$group <- dti_parameters_enrichment_asd$group <- 
  dti_parameters_enrichment_dd$group <- ica_nodes_enrichment_asd$group <- ica_nodes_enrichment_dd$group <- 'Significant'
regional_brain_volumes_enrichment_cond_asd$group <- regional_brain_volumes_enrichment_cond_dd$group <- dti_parameters_enrichment_cond_asd$group <- 
  dti_parameters_enrichment_cond_dd$group <- ica_nodes_enrichment_cond_asd$group <- ica_nodes_enrichment_cond_dd$group <- 'Conditionally significant'

# Put all results together
all_enrichment_result <- rbind.data.frame(
  regional_brain_volumes_enrichment_asd,
  dti_parameters_enrichment_asd,
  ica_nodes_enrichment_asd,
  regional_brain_volumes_enrichment_cond_asd,
  dti_parameters_enrichment_cond_asd,
  ica_nodes_enrichment_cond_asd,
  regional_brain_volumes_enrichment_dd,
  dti_parameters_enrichment_dd,
  ica_nodes_enrichment_dd,
  regional_brain_volumes_enrichment_cond_dd,
  dti_parameters_enrichment_cond_dd,
  ica_nodes_enrichment_cond_dd
)

# Apply the Haldane–Anscombe correction when the odds ratio is infinite, while retaining the original contingency-table counts.
all_enrichment_result$odds_ratio_adjusted <- all_enrichment_result$odds_ratio
idx <- is.infinite(all_enrichment_result$odds_ratio)
all_enrichment_result$odds_ratio_adjusted[idx] <- with(
  all_enrichment_result[idx, ],
  ((a + 0.5) * (d + 0.5)) /
    ((b + 0.5) * (c + 0.5))
)

# Adjust P values jointly across all trait-level enrichment tests
all_enrichment_result$p_fdr <- p.adjust(all_enrichment_result$p_hyper, method = "fdr")
all_enrichment_result$significant <- all_enrichment_result$p_fdr < 0.05

write.csv(all_enrichment_result, 'disease_rare_gene_enrichment_result.csv', row.names=F)