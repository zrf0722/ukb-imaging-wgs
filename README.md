# Whole-genome sequencing analyses of multimodal brain imaging phenotypes

This repository contains analysis code accompanying the manuscript **“Whole-genome sequencing analysis of multimodal brain imaging phenotypes.”** The study integrates whole-genome sequencing (WGS) and multimodal brain MRI data from approximately 80,000 UK Biobank participants to investigate the contributions of rare coding and noncoding regulatory variation to human brain structure and function.

The analyses cover 287 imaging-derived phenotypes (IDPs):

- 101 regional brain volumes;
- 110 diffusion tensor imaging (DTI) parameters; and
- 76 independent component analysis (ICA)-derived functional network traits.

The repository provides scripts for WGS preprocessing, heritability estimation, variant- and gene-based association testing, functional annotation enrichment analyses, and disease-gene enrichment analyses. 


## Repository contents

| File | Description |
| --- | --- |
| `WGS_preprocessing.sh` | Preprocesses and quality-controls WGS data for downstream analyses. |
| `WGS_heritability.sh` | Estimates WGS-based heritability and partitions heritability by allele frequency and functional annotation. |
| `summarize_WGS_heritability.R` | Summarizes and organizes WGS-based heritability results across imaging modalities and annotation categories. |
| `split_WGS_annotation_into_chunks.R` | Supports `WGS_heritability.sh` to divide chromosome-specific files into computational chunks. |
| `WGS_association_test_with_REGENIE.sh` | Performs variant-level and gene-based association tests using REGENIE. |
| `brain_functional_annotation_enrichment.R` | Tests for enrichment of significant variant-level associations in brain cell-type-specific ATAC-seq peaks. |
| `disease_rare_gene_enrichment.R` | Tests for enrichment of significant gene-based associations in disease-relevant gene sets. |
