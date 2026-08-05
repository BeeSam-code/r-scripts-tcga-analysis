# TCGA R Scripts - Analysis Guide

## Workflow Overview

This repository provides a modular workflow for TCGA data analysis:

```
Raw TCGA Data
    ↓
Data Processing
    ├── 01_rnaseq_data_processing.R
    ├── 02_cnv_data_processing.R
    └── 03_methylation_data_processing.R
    ↓
Data Integration & QC
    ↓
Analysis
    ├── 01_snp_ancestry_analysis.R
    └── 02_eqtl_identification.R
    ↓
Results & Visualization
```

## Data Processing Phase

### Step 1: Load Libraries

```r
source('utils/load_libraries.R')
```

This loads all required dependencies and handles installation if needed.

### Step 2: Process RNA-seq Data

```r
source('data-processing/01_rnaseq_data_processing.R')

# Main function
rnaseq_data <- load_and_process_rnaseq(
  data_dir = './data/rnaseq/',
  pattern = '.augmented_star_gene_counts.tsv$',
  min_expression = 0.1,
  min_prevalence = 0.8
)
```

**Outputs:**
- `COADREAD_rnaSeq_TPMonly_Jan2024.csv` - TPM normalized counts
- `COADREAD_tcga_InverseNormalTPMgenes_withShapirosig_Jan2024.csv` - INT-transformed expression

**QC Steps:**
- Removes genes with <50% non-zero values
- Filters genes by expression threshold
- Applies rank-based inverse normal transformation (INT)
- Shapiro-Wilk test for normality (p > 0.001)

### Step 3: Process CNV Data

```r
source('data-processing/02_cnv_data_processing.R')

cnv_data <- load_and_process_cnv(
  data_dir = './data/cnv/',
  pattern = '.v36.tsv$'
)
```

**Outputs:**
- `gene_level_CNV_COADREAD_Jan2024.csv` - Gene-level copy number

**Processing:**
- Reads gene-level CNV files
- Aggregates multiple samples
- Removes duplicate gene columns

### Step 4: Process Methylation Data

```r
source('data-processing/03_methylation_data_processing.R')

methyl_data <- load_and_process_methylation(
  data_dir = './data/methylation/',
  pattern = '.methylation_array.sesame.level3betas.txt$'
)
```

**Outputs:**
- `COADREADaverage_methylation_beta_values_perGene_Jan2024.csv` - Gene-level methylation

**Processing:**
- Reads methylation probe-level data
- Aggregates to gene level
- Filters probes with >40% missing data
- Normalizes sample IDs

## Analysis Phase

### Step 5: SNP Ancestry Analysis

```r
source('analysis/01_snp_ancestry_analysis.R')

# VCF to GDS conversion
snpgdsVCF2GDS(
  vcf.fn = './data/vcf/merged_COADREAD_1KGP.vcf.gz',
  gds.fn = 'merged_tcgaCOADREAD_1kgp_Jan2024.gds'
)

# PCA for ancestry inference
ancestry_results <- perform_ancestry_pca(
  gds_file = 'merged_tcgaCOADREAD_1kgp_Jan2024.gds',
  sample_metadata = clinical_data
)
```

**Outputs:**
- PCA plot: `SNPrelate_COADREAD_PCAJan2024.png`
- Ancestry estimates table
- Admixture proportions

**Analysis:**
- Converts VCF to efficient GDS format
- Performs PCA with reference populations (CEU, YRI, JPT, CHB)
- Estimates ancestry proportions
- Classifies samples as European, African, or admixed

### Step 6: eQTL Identification

```r
source('analysis/02_eqtl_identification.R')

# Single gene eQTL analysis
eqtl_results <- identify_eqtls_single_gene(
  gene_name = 'POLB',
  genotype_data = genotype_matrix,
  expression_data = rnaseq_data,
  methylation_data = methyl_data,
  cnv_data = cnv_data,
  covariates = clinical_covariates,
  cis_distance = 1e6,
  pvalue_threshold = 0.01
)

# All genes eQTL analysis
eqtl_all <- identify_eqtls_all_genes(
  gene_list = common_genes,
  genotype_data = genotype_matrix,
  expression_data = rnaseq_data,
  covariates = clinical_covariates,
  output_dir = './results/eqtl/'
)
```

**Outputs:**
- `[gene]_eQTLnov23_cis.txt` - cis-eQTL results
- `[gene]_eQTLnov23_trans.txt` - trans-eQTL results

**Analysis:**
- Uses MatrixEQTL for efficient large-scale testing
- Includes covariates: age, gender, genetic ancestry (PCs), methylation, CNV
- Separate cis/trans analysis with distance threshold
- Multiple testing correction via permutation

## Data Integration

### Matching Samples Across Datasets

```r
# Find common samples across all data types
common_samples <- intersect_samples(
  rnaseq_samples = colnames(rnaseq_data),
  cnv_samples = colnames(cnv_data),
  methylation_samples = colnames(methyl_data),
  genotype_samples = colnames(genotype_matrix),
  ancestry_samples = colnames(ancestry_pcs)
)

# Filter all datasets to common samples
rnaseq_final <- rnaseq_data[, common_samples]
cnv_final <- cnv_data[, common_samples]
methyl_final <- methyl_data[, common_samples]
genotype_final <- genotype_matrix[, common_samples]
```

## Output Organization

```
results/
├── processed_data/
│   ├── rnaseq_normalized.csv
│   ├── cnv_gene_level.csv
│   ├── methylation_gene_level.csv
│   └── sample_metadata.csv
├── ancestry/
│   ├── pca_plot.png
│   ├── ancestry_estimates.csv
│   └── admixture_plot.png
└── eqtl/
    ├── gene1_cis_eqtls.txt
    ├── gene1_trans_eqtls.txt
    └── eqtl_summary.txt
```

## Troubleshooting

### Common Issues

**Issue**: Memory error with large datasets
- **Solution**: Process by chromosome or use data.table chunking

**Issue**: Sample ID mismatches across datasets
- **Solution**: Use `standardize_sample_ids()` function with mapping file

**Issue**: Missing values in covariates
- **Solution**: Use `impute_covariates()` or filter affected samples

### Performance Tips

1. Process data in batches if dataset is >10GB
2. Use data.table for faster operations on large files
3. Pre-filter genes to reduce computational load
4. Run eQTL analysis in parallel by chromosome

## References

- MatrixEQTL: https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
- SNPRelate: https://bioconductor.org/packages/release/bioc/html/SNPRelate.html
- DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
