#!/usr/bin/env Rscript

# Representative code for dividing a chromosome-specific SNP information file
# into computational chunks for PLINK extraction and MPH GRM construction.
#
# This script accompanies:
#   prepare_WGS_annotations_and_GRMs.sh

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

chromosome <- args[1]
snp_info_file <- args[2]

snp_info <- fread(snp_info_file)

chunk_size <- 800000L
n_variants <- nrow(snp_info)
chunk_starts <- seq(1L, n_variants, by = chunk_size)

for (part in seq_along(chunk_starts)) {

  start_row <- chunk_starts[part]
  end_row <- min(start_row + chunk_size - 1L, n_variants)

  snp_info_chunk <- snp_info[start_row:end_row]

  # Save the annotation file used by MPH.
  fwrite(snp_info_chunk, file = sprintf("snp_info_chr%s_part%d_coding_noncoding.csv", chromosome, part), quote = FALSE)

  # Save the variant list used for PLINK extraction.
  fwrite(snp_info_chunk[, .(SNP)], file = sprintf("chr%s_part%d_all.snps", chromosome, part), col.names = FALSE, quote = FALSE)
}


