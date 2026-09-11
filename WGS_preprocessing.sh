#!/usr/bin/env bash

# Representative UKB-RAP workflow for preprocessing WGS pVCFs.
#
# Main steps:
#   1. Trim chromosome-specific pVCF shards using an AAScore threshold.
#   2. Normalize trimmed VCFs against the GRCh38 reference.
#   3. Concatenate ordered normalized VCF shards by chromosome.
#   4. Convert chromosome-level VCFs to BGEN v1.2.
#
# Project-specific paths have been replaced with placeholders. This script is
# representative code and is not intended to serve as a directly executable
# end-to-end pipeline.

set -euo pipefail


# -----------------------------------------------------------------------------
# Project paths
# -----------------------------------------------------------------------------

PROJECT_ROOT="/PATH_TO_PROJECT"

RAW_LIST_DIR="/PATH_TO_RAW_PVCF_LISTS"
TRIMMED_LIST_DIR="/PATH_TO_TRIMMED_VCF_LISTS"
NORMALIZED_LIST_DIR="/PATH_TO_NORMALIZED_VCF_LISTS"

TRIMMED_VCF_DIR="${PROJECT_ROOT}/trimmed_VCF"
NORMALIZED_VCF_DIR="${PROJECT_ROOT}/normalized_VCF"
MERGED_VCF_DIR="${PROJECT_ROOT}/merged_VCF"
BGEN_DIR="${PROJECT_ROOT}/BGEN"

REFERENCE_FASTA="/PATH_TO/GRCh38_full_analysis_set_plus_decoy_hla.fa"

TRIMMER_APPLET="/PATH_TO/vcf_trimmer"

CHROMOSOMES=($(seq 1 22))


# -----------------------------------------------------------------------------
# Step 1. Trim pVCF shards
# -----------------------------------------------------------------------------
# Each input list contains the raw pVCF shards for one chromosome.
# In the original workflow, each chromosome list was divided into multiple
# batches and submitted as parallel UKB-RAP jobs.
#
# Variants with INFO/AAScore <= 0.5 were removed, and FORMAT fields not needed
# for downstream analyses were excluded.

for chr in "${CHROMOSOMES[@]}"; do
  dx run "${TRIMMER_APPLET}" \
    -ivcf_file_list="${RAW_LIST_DIR}/chr${chr}_vcf_list" \
    -ifile_label="trimmed" \
    -ioutput_dir="${TRIMMED_VCF_DIR}/chr${chr}" \
    -iqc_thresholds="INFO/AAScore>0.5" \
    -ifields_to_remove="FORMAT/FT,FORMAT/AD,FORMAT/MD,FORMAT/DP,FORMAT/RA,FORMAT/PP,FORMAT/GQ,FORMAT/PL" \
    --tag="Trim WGS pVCFs chr${chr}" \
    --yes
done


# -----------------------------------------------------------------------------
# Step 2. Normalize trimmed VCF shards
# -----------------------------------------------------------------------------
# Each chromosome-specific text file lists the trimmed VCF shards to normalize.
# Normalization was performed against the GRCh38 reference using bcftools norm.

for chr in "${CHROMOSOMES[@]}"; do
  while IFS= read -r input_vcf; do
    [[ -z "${input_vcf}" ]] && continue

    input_name="$(basename "${input_vcf}")"
    output_name="${input_name%.vcf.gz}_normalized.vcf.gz"

    bcftools norm \
      "${input_vcf}" \
      -f "${REFERENCE_FASTA}" \
      -Oz \
      -o "${NORMALIZED_VCF_DIR}/chr${chr}/${output_name}"

  done < "${TRIMMED_LIST_DIR}/chr${chr}_vcf_list"
done


# -----------------------------------------------------------------------------
# Step 3. Concatenate normalized VCF shards by chromosome
# -----------------------------------------------------------------------------
# Each normalized-VCF list must contain nonoverlapping chromosome shards in
# genomic order, with compatible headers and identical sample ordering.

for chr in "${CHROMOSOMES[@]}"; do
  bcftools concat \
    -f "${NORMALIZED_LIST_DIR}/chr${chr}_vcf_list" \
    -Oz \
    -o "${MERGED_VCF_DIR}/chr${chr}_merged_normalized.vcf.gz"
done


# -----------------------------------------------------------------------------
# Step 4. Convert chromosome-level VCFs to BGEN v1.2
# -----------------------------------------------------------------------------

for chr in "${CHROMOSOMES[@]}"; do
  plink2 \
    --vcf "${MERGED_VCF_DIR}/chr${chr}_merged_normalized.vcf.gz" \
    --export bgen-1.2 bits=8 ref-first \
    --out "${BGEN_DIR}/chr${chr}_merged_normalized"
done