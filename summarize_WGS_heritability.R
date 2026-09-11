#!/usr/bin/env Rscript

# Representative code for summarizing WGS-based heritability estimates
# partitioned by coding status and allele-frequency category.
#
# Project-specific input paths and variant counts have been replaced with
# placeholders.
#
# This script is not intended to serve as a standalone executable pipeline.

library(data.table)

phenotype_file <- "PATH_TO_DTI_PHENOTYPE_FILE"

phenotypes <- fread(phenotype_file)
traits <- colnames(phenotypes)[-1]

# Numbers of coding-rare, coding-common, noncoding-rare, and noncoding-common variants used in the enrichment calculations
M_coding_rare <- Number_of_coding_rare_variants
M_coding_common <- Number_of_coding_common_variants
M_noncoding_rare <- Number_of_noncoding_rare_variants
M_noncoding_common <- Number_of_noncoding_common_variants

M_common <- M_coding_common + M_noncoding_common 
M_rare <- M_coding_rare + M_noncoding_rare 
M_noncoding <- M_noncoding_common + M_noncoding_rare 
M_coding <- M_coding_common + M_coding_rare 
M_total <- M_coding_common + M_noncoding_common + M_coding_rare + M_noncoding_rare

# get se of coding common/rare enrichment
get_ratio_enrichment_se <- function(x, y, var_x, var_y, cov_xy, c = 1) {
  var_e <- (c / y)^2 * var_x + (c * x / y^2)^2 * var_y - 2 * (c / y) * (c * x / y^2) * cov_xy
  sqrt(var_e)
}

cc  <- "PATH_TO_CODING_COMMON_GRM_PREFIX"
cr  <- "PATH_TO_CODING_RARE_GRM_PREFIX"
ncc <- "PATH_TO_NONCODING_COMMON_GRM_PREFIX"
ncr <- "PATH_TO_NONCODING_RARE_GRM_PREFIX"

h2 <- matrix(NA, ncol=4, nrow=0)
h2_se <- matrix(NA, ncol=1, nrow=0)
h2_se_rare <- matrix(NA, ncol=1, nrow=0)
h2_se_common <- matrix(NA, ncol=1, nrow=0)
enrichment_se <- matrix(NA, ncol=3, nrow=0)
for (trait in traits) {
  res <- fread(paste0('heritability_reml_', trait, '.mq.vc.csv'))

  h2 <- rbind.data.frame(h2, res$pve[1:4])
  h2_se <- rbind.data.frame(h2_se, sqrt(sum(res$seP[1:4]^2)))
  h2_se_rare <- rbind.data.frame(h2_se_rare, sqrt(sum(res$seP[c(2,4)]^2)))
  h2_se_common <- rbind.data.frame(h2_se_common, sqrt(sum(res$seP[c(1,3)]^2)))

  vc_names <- res$vc_name
  n_vc <- length(vc_names)
  start_pve_cov <- ncol(res) - n_vc + 1
  S <- as.matrix(res[, start_pve_cov:ncol(res), with = FALSE])
  rownames(S) <- vc_names
  colnames(S) <- vc_names
  se_enrichment_coding <- get_ratio_enrichment_se(x = sum(res$pve[c(1, 2)]), y = sum(res$pve[c(1, 2, 3, 4)]),
    var_x = S[cc, cc] + S[cr, cr] + 2 * S[cc, cr], var_y = sum(S[c(cc, cr, ncc, ncr), c(cc, cr, ncc, ncr)]),
    cov_xy = sum(S[c(cc, cr), c(cc, cr, ncc, ncr)]), c = M_total / (M_coding_common + M_coding_rare)
  )
  se_enrichment_coding_common <- get_ratio_enrichment_se(x = res$pve[1], y = sum(res$pve[c(1,3)]),
    var_x = S[cc, cc], var_y = S[cc, cc] + S[ncc, ncc] + 2 * S[cc, ncc],
    cov_xy = S[cc, cc] + S[cc, ncc], c = M_common / M_coding_common
  )
  se_enrichment_coding_rare <- get_ratio_enrichment_se(x = res$pve[2], y = sum(res$pve[c(2,4)]),
    var_x = S[cr, cr], var_y = S[cr, cr] + S[ncr, ncr] + 2 * S[cr, ncr],
    cov_xy = S[cr, cr] + S[cr, ncr], c = M_rare / M_coding_rare
  )
  enrichment_se <- rbind.data.frame(enrichment_se, c(se_enrichment_coding, se_enrichment_coding_rare, se_enrichment_coding_common))
}
h2 <- cbind.data.frame(traits, h2, h2_se, h2_se_rare, h2_se_common)
colnames(h2) <- c('trait', 'h2_coding_common', 'h2_coding_rare', 'h2_noncoding_common', 'h2_noncoding_rare', 'se_total', 'se_rare', 'se_common')
h2$h2_total <- h2$h2_coding_common + h2$h2_coding_rare + h2$h2_noncoding_common + h2$h2_noncoding_rare 
h2$h2_coding <- h2$h2_coding_common + h2$h2_coding_rare
h2$h2_noncoding <- h2$h2_noncoding_common + h2$h2_noncoding_rare
h2$h2_rare <- h2$h2_coding_rare + h2$h2_noncoding_rare
h2$h2_common <- h2$h2_coding_common + h2$h2_noncoding_common
h2$coding_frac <- h2$h2_coding / h2$h2_total

h2$enrichment_coding <- (h2$h2_coding / M_coding) / (h2$h2_total / M_total)
h2$enrichment_coding_common <- (h2$h2_coding_common / M_coding_common) / (h2$h2_common / M_common)
h2$enrichment_coding_rare   <- (h2$h2_coding_rare / M_coding_rare) / (h2$h2_rare / M_rare)

enrichment_se <- cbind.data.frame(traits, enrichment_se)
colnames(enrichment_se) <- c('trait', 'se_coding_enrichment', 'se_coding_rare_enrichment', 'se_coding_common_enrichment')
h2_all <- merge(h2, enrichment_se, by='trait', all.x=T)
h2_all$p_coding_enrichment <- 2 * pnorm(-abs((h2_all$enrichment_coding - 1) / h2_all$se_coding_enrichment))
h2_all$p_coding_rare_enrichment <- 2 * pnorm(-abs((h2_all$enrichment_coding_rare - 1) / h2_all$se_coding_rare_enrichment))
h2_all$p_coding_common_enrichment <- 2 * pnorm(-abs((h2_all$enrichment_coding_common - 1) / h2_all$se_coding_common_enrichment))

# Total h2 Wald test
h2_all$z_total <- h2_all$h2_total / h2_all$se_total
h2_all$p_total <- 2 * pnorm(-abs(h2_all$z_total))

# Rare h2 Wald test
h2_all$z_rare <- h2_all$h2_rare / h2_all$se_rare
h2_all$p_rare <- 2 * pnorm(-abs(h2_all$z_rare))

# Significance required a Bonferroni-corrected total heritability test across all imaging traits
# and nominal significance for rare-variant heritability.
n_total_traits <- 110 + 101 + 76
h2_all$sig <- h2_all$p_total < 0.05 / n_total_traits & h2_all$p_rare < 0.05
write.csv(h2_all, 'heritability_wgs.csv', row.names =  FALSE)
