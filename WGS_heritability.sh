#!/usr/bin/env bash

# Representative commands for annotating WGS variants and constructing
# annotation-stratified genetic relationship matrices (GRMs) for the
# WGS-based heritability analyses described in the manuscript.
#
# Variants were classified by coding status and minor allele frequency (MAF) into:
#   1. coding common
#   2. coding rare
#   3. noncoding common
#   4. noncoding rare
#
# Project-specific paths, chromosome indices, chunk indices, and computing-
# environment settings are represented by shell variables.
#
# This script is not intended to serve as a standalone executable pipeline.


# ==============================================================================
# Variables used in the representative commands
# ==============================================================================

# CHR:
#   Chromosome number, from 1 to 22.
#
# PART:
#   Chromosome chunk number.
#
# ANNOTATION:
#   One of:
#     coding_common
#     coding_rare
#     noncoding_common
#     noncoding_rare
#
# VEP_CACHE_DIR:
#   Path to the local VEP cache.
#
# WGS_BFILE_PREFIX:
#   Prefix of the chromosome-specific QC-filtered PLINK dataset.
#
# PHENOTYPE_FILE:
#   Phenotype file.
#
# COVARIATE_FILE:
#   Covariate file.
# 
# UNRELATED_SAMPLE_FILE:
#   Sample IDs included in the unrelated analysis set.
#
# GRM_LIST:
#   File listing the four unrelated annotation-specific GRM prefixes.

VEP_CACHE_DIR="PATH_TO_VEP_CACHE"
WGS_BFILE_PREFIX="PATH_TO_CHROMOSOME_SPECIFIC_WGS_PLINK_PREFIX"
WGS_VCF="PATH_TO_CHROMOSOME_SPECIFIC_WGS_VCF"

PHENOTYPE_FILE="PATH_TO_PHENOTYPE_FILE"
COVARIATE_FILE="PATH_TO_COVARIATE_FILE"
UNRELATED_SAMPLE_FILE="PATH_TO_UNRELATED_SAMPLE_LIST"
GRM_LIST="PATH_TO_REML_GRM_LIST"


# ==============================================================================
# 1. Annotate WGS variants using VEP
# ==============================================================================

# VEP was run using the offline cache and one selected consequence per variant.
# The downstream parsing assumes that the VEP CSQ format begins with:
#
#   Allele|Consequence|IMPACT|...
#
# Coding variants included synonymous and protein-altering coding consequences.

vep \
  -i "${WGS_VCF}"
  -o "chr${CHR}.vep.vcf" \
  --format vcf \
  --vcf \
  --pick \
  --offline \
  --cache \
  --dir_cache "${VEP_CACHE_DIR}" \
  --fork 4 \
  --force


# Compress and index the VEP-annotated VCF.
bgzip -c "chr${CHR}.vep.vcf" > chr${CHR}.vep.vcf.gz
tabix -p vcf chr${CHR}.vep.vcf.gz


# ==============================================================================
# 2. Calculate MAF and classify variants as coding or noncoding
# ==============================================================================

# Calculate alternate-allele frequencies.
plink2 \
  --bfile "${WGS_BFILE_PREFIX}" \
  --freq cols=+altfreq \
  --out "chr${CHR}.wgs_frequency"


# Create a two-column file containing variant ID and MAF for variants with MAF > 0.0001.
awk -v output_file="chr${CHR}.mafgt0p0001.tsv" '
NR == 1 {next}
{
    af = $6; maf = (af < 0.5 ? af : 1 - af);
    if (maf > 0.0001) {
      print $2 "\t" maf > output_file
    }
}
' "chr${CHR}.wgs_frequency.afreq"


# Extract coding status from the selected VEP consequence.
bcftools query -f '%CHROM:%POS:%REF:%ALT\t%INFO/CSQ\n' "chr${CHR}.vep.vcf.gz" |
awk -F '\t' '
BEGIN {
    split(
        "missense_variant,stop_gained,stop_lost,frameshift_variant," \
        "splice_acceptor_variant,splice_donor_variant,start_lost," \
        "inframe_insertion,inframe_deletion,synonymous_variant",
        coding_terms,
        ","
    )
    for (i in coding_terms) {
        coding[coding_terms[i]] = 1
    }
    OFS = "\t"
}

{
    variant_id = $1
    csq = $2

    # Consequence is assumed to be the second field of the VEP CSQ annotation:
    # Allele|Consequence|IMPACT|...
    split(csq, csq_fields, "\\|")
    consequence = csq_fields[2]

    # A selected consequence may contain multiple terms joined by "&".
    split(consequence, consequence_terms, "&")

    is_coding = 0

    for (i in consequence_terms) {
        if (consequence_terms[i] in coding) {
            is_coding = 1
            break
        }
    }

    print variant_id, is_coding
}
' > "chr${CHR}.coding_status.tsv"


# Combine coding status and MAF into the four annotation categories.
awk 'BEGIN{OFS=","; print "SNP,coding_rare,coding_common,noncoding_rare,noncoding_common"}

NR == FNR{maf[$1]=$2; next}
{
    variant_id = $1
    is_coding = $2

    if (!(variant_id in maf)) next;

    variant_maf = maf[variant_id]

    coding_rare = (is_coding == 1 && variant_maf > 0.0001 && variant_maf <= 0.01) ? 1 : ""
    coding_common = (is_coding == 1 && variant_maf > 0.01) ? 1 : ""
    noncoding_rare = (is_coding == 0 && variant_maf > 0.0001 && variant_maf <= 0.01) ? 1 : ""
    noncoding_common = (is_coding == 0 && variant_maf > 0.01) ? 1 : ""
    
    print variant_id, coding_rare, coding_common, noncoding_rare, noncoding_common
}' "chr${CHR}.mafgt0p0001.tsv" "chr${CHR}.coding_status.tsv" > "chr${CHR}.snp_info.coding.noncoding.csv"

# Remove duplicate complete rows while retaining the header.
head -1 chr${CHR}.snp_info.coding.noncoding.csv > chr${CHR}.snp_info.coding.noncoding.cleaned.csv
tail -n +2 chr${CHR}.snp_info.coding.noncoding.csv | awk '!seen[$0]++' >> chr${CHR}.snp_info.coding.noncoding.cleaned.csv


# ==============================================================================
# 4. Divide each chromosome into computational chunks
# ==============================================================================

# The companion R script generates:
#
#   snp_info_chrN_partM_coding_noncoding.csv
#   chrN_partM_all.snps
#
# Each chunk contains up to 800,000 variants.

Rscript split_WGS_annotation_into_chunks.R ${CHR} chr${CHR}.snp_info.coding_noncoding.cleaned.csv


# ==============================================================================
# 5. Extract PLINK data for each chromosome chunk
# ==============================================================================

plink \
  --bfile "${WGS_BFILE_PREFIX}" \
  --extract "chr${CHR}_part${PART}_all.snps" \
  --make-bed \
  --out "chr${CHR}_part${PART}_wgs"


# ==============================================================================
# 6. Construct annotation-specific GRMs for each chunk
# ==============================================================================

# This command was run separately for each annotation category:
#
#   coding_common
#   coding_rare
#   noncoding_common
#   noncoding_rare

mph \
  --bfile "chr${CHR}_part${PART}_wgs" \
  --snp_info_file "snp_info_chr${CHR}_part${PART}_coding_noncoding.csv" \
  --make_grm \
  --snp_weight_name "${ANNOTATION}" \
  --num_threads 10 \
  --output_file "grm_chr${CHR}_part${PART}_${ANNOTATION}"


# ==============================================================================
# 7. Merge chunk-level GRMs within each chromosome
# ==============================================================================

# The GRM list contains the chunk-level GRM prefixes for one chromosome and
# one annotation category.

mph \
  --grm_list "grm_list_chr${CHR}_${ANNOTATION}.txt" \
  --merge_grms \
  --num_threads 10 \
  --output_file "grm_chr${CHR}_${ANNOTATION}"


# ==============================================================================
# 8. Merge chromosome-level GRMs across chromosomes 1-22
# ==============================================================================

# The genome-wide GRM list contains the chromosome-level GRM prefixes for
# one annotation category.

mph \
  --grm_list "grm_list_chr1_22_${ANNOTATION}.txt" \
  --merge_grms \
  --num_threads 10 \
  --output_file "grm_allchr_${ANNOTATION}"


# ==============================================================================
# 9. Subset each annotation-specific GRM to unrelated participants
# ==============================================================================

# This command was run separately for each annotation category.

mph \
  --subset_grm \
  --binary_grm_file "grm_allchr_${ANNOTATION}" \
  --keep "${UNRELATED_SAMPLE_FILE}" \
  --num_threads 10 \
  --output_file "grm_allchr_${ANNOTATION}_unrelated"


# ==============================================================================
# 10. Estimate annotation-partitioned heritability using REML
# ==============================================================================

# The GRM list supplied to MPH contained the following genome-wide,
# unrelated-sample GRM prefixes:
#
#   grm_allchr_coding_common_unrelated
#   grm_allchr_coding_rare_unrelated
#   grm_allchr_noncoding_common_unrelated
#   grm_allchr_noncoding_rare_unrelated

mph \
  --reml \
  --save_memory \
  --grm_list "${GRM_LIST}" \
  --phenotype_file "${PHENOTYPE_FILE}" \
  --trait "${TRAIT}" \
  --covariate_file "${COVARIATE_FILE}" \
  --num_threads 10 \
  --output_file "heritability_reml_${TRAIT}"


  