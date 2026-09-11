#!/usr/bin/env bash

# Representative commands for the variant-level and gene-based association
# analyses described in the manuscript "Whole-genome sequencing analysis of 
# multimodal brain imaging phenotypes".
#
# These commands document the principal software options and analysis settings.
# Project-specific file paths, phenotype names, chromosome indices, and
# computing-environment settings are represented by shell variables and
# must be adapted before execution.
#
# This script is not intended to serve as a standalone executable pipeline.


# ==============================================================================
# Genotype data quality control
# ==============================================================================

# Genotyping-array QC for REGENIE Step 1
plink \
  --bfile "${array_genotype_prefix}" \
  --keep "${sample_inclusion_file}" \
  --mind 0.1 \
  --geno 0.1 \
  --hwe 1e-15 \
  --maf 0.01 \
  --mac 20 \
  --make-bed \
  --out "${array_qc_prefix}"


# WGS QC for REGENIE Step 2
plink2 \
  --autosome \
  --bgen "${chrN_bgen}" ref-first \
  --sample "${chrN_sample}" \
  --keep "${sample_inclusion_file}" \
  --geno 0.1 \
  --hwe 1e-30 \
  --export bgen-1.2 bits=8 ref-first \
  --rm-dup force-first list \
  --out "${chrN_wgs_qc_prefix}"


# ==============================================================================
# Association analyses with REGENIE
# ==============================================================================

# REGENIE Step 1
regenie \
  --step 1 \
  --bed "${array_qc_prefix}" \
  --phenoFile "${phenotype_file}" \
  --covarFile "${covariate_file}" \
  --lowmem \
  --lowmem-prefix "${temporary_prefix}" \
  --apply-rint \
  --bsize 1000 \
  --out "${step1_output_prefix}"


# REGENIE Step 2: variant-level association testing
regenie \
  --step 2 \
  --bgen "${chrN_wgs_qc_prefix}.bgen" \
  --sample "${chrN_wgs_qc_prefix}.sample" \
  --bgi "${chrN_wgs_qc_prefix}.bgen.bgi" \
  --ref-first \
  --phenoFile "${phenotype_file}" \
  --covarFile "${covariate_file}" \
  --phenoCol "${phenotype_name}" \
  --pred "${step1_pred_list}" \
  --apply-rint \
  --minMAC 1 \
  --bsize 400 \
  --gz \
  --out "${variant_test_output_prefix}"
  

# REGENIE Step 2: gene-based association testing
regenie \
  --step 2 \
  --bgen "${chrN_wgs_qc_prefix}.bgen" \
  --sample "${chrN_wgs_qc_prefix}.sample" \
  --bgi "${chrN_wgs_qc_prefix}.bgen.bgi" \
  --ref-first \
  --phenoFile "${phenotype_file}" \
  --covarFile "${covariate_file}" \
  --phenoCol "${phenotype_name}" \
  --pred "${step1_pred_list}" \
  --apply-rint \
  --minMAC 1 \
  --anno-file "${annotation_file}" \
  --set-list "${set_list_file}" \
  --mask-def "${mask_definition_file}" \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --aaf-bins 0.01 \
  --bsize 400 \
  --gz \
  --out "${gene_test_output_prefix}"


# ==============================================================================
# Conditional analyses
# ==============================================================================

# Convert GWAS summary statistics to GCTA-COJO format
#
# The selected input columns correspond to:
# SNP, effect allele, other allele, allele frequency, effect estimate,
# standard error, P value, and sample size.
awk 'BEGIN{OFS="\t"}
     NR==1 {print "SNP","A1","A2","freq","b","se","p","N"; next}
     {print $2,$4,$5,$7,$8,$9,$10,$6}' \
  "${gwas_summary_statistics}" > "${cojo_input_file}"


# Identify independent significant GWAS variants using GCTA-COJO
gcta64 \
  --bfile "${ld_reference_prefix}" \
  --cojo-file "${cojo_input_file}" \
  --chr "${chromosome}" \
  --cojo-slct \
  --cojo-collinear 0.9 \
  --cojo-wind 10000 \
  --out "${cojo_output_prefix}"


# REGENIE Step 2: variant-level association testing conditioned on
# independently associated GWAS variants
regenie \
  --step 2 \
  --bgen "${chrN_wgs_qc_prefix}.bgen" \
  --sample "${chrN_wgs_qc_prefix}.sample" \
  --bgi "${chrN_wgs_qc_prefix}.bgen.bgi" \
  --ref-first \
  --phenoFile "${phenotype_file}" \
  --covarFile "${covariate_file}" \
  --phenoCol "${phenotype_name}" \
  --pred "${step1_pred_list}" \
  --apply-rint \
  --minMAC 1 \
  --bsize 400 \
  --gz \
  --extract "${variants_for_conditional_analysis}" \
  --condition-list "${independent_gwas_variant_list}" \
  --condition-file "bed,${conditioning_genotype_prefix}" \
  --out "${conditional_variant_test_output_prefix}"


# REGENIE Step 2: gene-based association testing conditioned on
# independently associated GWAS variants
regenie \
  --step 2 \
  --bgen "${chrN_wgs_qc_prefix}.bgen" \
  --sample "${chrN_wgs_qc_prefix}.sample" \
  --bgi "${chrN_wgs_qc_prefix}.bgen.bgi" \
  --ref-first \
  --phenoFile "${phenotype_file}" \
  --covarFile "${covariate_file}" \
  --phenoCol "${phenotype_name}" \
  --pred "${step1_pred_list}" \
  --apply-rint \
  --minMAC 1 \
  --anno-file "${annotation_file}" \
  --set-list "${set_list_file}" \
  --mask-def "${mask_definition_file}" \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --aaf-bins 0.01 \
  --bsize 400 \
  --gz \
  --extract-sets "${genes_for_conditional_analysis}" \
  --condition-list "${independent_gwas_variant_list}" \
  --condition-file "bed,${conditioning_genotype_prefix}" \
  --out "${conditional_gene_test_output_prefix}"