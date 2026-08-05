# TCGA R Scripts - Script Descriptions

## Data Processing Scripts

### 01_rnaseq_data_processing.R

**Purpose**: Load, process, and normalize RNA-seq count data from TCGA

**Main Functions**:
- `load_raw_rnaseq_files()` - Reads multiple RNA-seq count files
- `filter_genes_by_expression()` - Removes lowly expressed genes
- `apply_normalization()` - TPM normalization and rank-based inverse normal transformation
- `qc_expression_data()` - Quality control and normality testing

**Input**: TCGA RNA-seq augmented STAR gene count files (.tsv)

**Output**: 
- TPM-normalized counts
- INT-transformed expression data
- QC metrics

**Key Processing Steps**:
1. Read count files with gene IDs and sample metadata
2. Remove genes with >50% zero counts
3. Filter genes with low expression (default: <0.1 TPM in 80% of samples)
4. Apply rank-based inverse normal transformation
5. Test normality with Shapiro-Wilk test
6. Retain genes passing normality threshold (p > 0.001)

**Parameters**:
- `data_dir`: Directory containing count files
- `pattern`: File pattern to match (e.g., ".augmented_star_gene_counts.tsv$")
- `min_expression`: Minimum expression threshold
- `min_prevalence`: Minimum proportion of samples with detectable expression
- `normality_threshold`: p-value threshold for Shapiro-Wilk test

---

### 02_cnv_data_processing.R

**Purpose**: Load and aggregate gene-level copy number variation data

**Main Functions**:
- `load_cnv_files()` - Reads gene-level CNV files from GDC
- `aggregate_cnv_data()` - Combines multiple sample CNV files
- `remove_duplicate_columns()` - Eliminates redundant gene name columns
- `standardize_sample_ids()` - Normalizes TCGA sample barcodes

**Input**: TCGA gene-level CNV files (.v36.tsv)

**Output**: 
- Gene-level CNV matrix (genes × samples)
- Standardized sample naming

**Key Processing Steps**:
1. Read gene-level copy number files
2. Extract sample ID from file path using UUIDtoBarcode conversion
3. Combine all samples using cbindX for efficient column binding
4. Remove redundant gene name columns (keep only first)
5. Standardize sample IDs to TCGA format (TCGA-XX-XXXX)
6. Export as CSV with gene names as row names

**Parameters**:
- `data_dir`: Directory containing CNV files
- `pattern`: File pattern to match (e.g., ".v36.tsv$")
- `min_samples`: Minimum samples with non-missing data

---

### 03_methylation_data_processing.R

**Purpose**: Process methylation array data and aggregate to gene level

**Main Functions**:
- `load_methylation_files()` - Reads probe-level methylation data
- `annotate_probes()` - Maps probes to genes using HM450 manifest
- `aggregate_to_gene_level()` - Averages methylation across probes per gene
- `filter_missing_data()` - Removes probes/genes with excessive missing values
- `normalize_beta_values()` - Standardizes beta-value distributions

**Input**: 
- TCGA methylation array files (.sesame.level3betas.txt)
- HM450 reference manifest (hg38)

**Output**: 
- Gene-level methylation matrix (genes × samples)
- Beta-value normalized data

**Key Processing Steps**:
1. Read probe-level methylation beta values
2. Annotate probes using HM450 manifest to map to genes
3. Handle genes with multiple probe sets (average across probes)
4. Filter probes with >40% missing data
5. Filter genes with >40% missing data across samples
6. Normalize sample IDs to TCGA format
7. Prefix gene names with 'methyl_' to avoid conflicts with other data types

**Parameters**:
- `data_dir`: Directory containing methylation files
- `pattern`: File pattern to match
- `manifest_file`: Path to HM450 annotation manifest
- `max_missing_rate`: Maximum proportion of missing data (default: 0.4)

---

## Analysis Scripts

### 01_snp_ancestry_analysis.R

**Purpose**: Perform SNP-based ancestry inference using PCA and admixture estimation

**Main Functions**:
- `convert_vcf_to_gds()` - Converts VCF to efficient GDS format
- `perform_pca_analysis()` - Performs PCA on genotype data
- `estimate_ancestry()` - Calculates ancestry proportions
- `plot_ancestry_results()` - Creates PCA and admixture plots

**Input**: 
- VCF file with TCGA and 1000 Genomes samples merged
- Clinical data with self-reported race/ethnicity

**Output**: 
- PCA plot (PNG)
- Ancestry classification table
- Admixture proportion estimates

**Key Analysis Steps**:
1. Convert VCF to GDS format for efficient access
2. Sample 1000 SNPs for initial PCA
3. Perform PCA across all samples
4. Extract principal components (PC1-PC3)
5. Map reference populations (CEU, YRI, JPT, CHB) to PCA space
6. Classify TCGA samples based on PC positions
7. Estimate admixture proportions using SNPRelate
8. Generate publication-quality plots

**Parameters**:
- `vcf_file`: Path to merged VCF file
- `gds_file`: Output GDS file path
- `reference_pops`: List of reference populations
- `num_snps`: Number of SNPs to sample (default: 1000)
- `ancestry_threshold`: Threshold for ancestry classification (default: 0.60)

---

### 02_eqtl_identification.R

**Purpose**: Identify expression quantitative trait loci (eQTLs) using MatrixEQTL

**Main Functions**:
- `prepare_eqtl_matrices()` - Formats genotype, expression, and covariate data
- `identify_eqtls_single_gene()` - Tests single gene for cis/trans-eQTLs
- `identify_eqtls_all_genes()` - Performs genome-wide eQTL analysis
- `annotate_eqtl_results()` - Adds genomic annotations to results
- `plot_eqtl_manhattan()` - Creates Manhattan plots

**Input**: 
- Genotype matrix (SNPs × samples, biallelic coding)
- Expression matrix (genes × samples, normalized)
- Covariate matrix (age, gender, genetic PCs, methylation, CNV)
- SNP and gene position files

**Output**: 
- cis-eQTL results file
- trans-eQTL results file
- Visualization plots

**Key Analysis Steps**:
1. Format data using SlicedData objects for memory efficiency
2. Load SNP positions and gene positions
3. Define covariates (age, gender, genetic ancestry PCs, methylation, CNV)
4. Specify cis-distance threshold (default: 1Mb)
5. Test for association between SNPs and gene expression
6. Perform separate cis (SNPs near gene) and trans (genome-wide) analysis
7. Apply p-value thresholds and FDR correction
8. Output significant associations with effect sizes and p-values

**Parameters**:
- `cis_distance`: Maximum distance for cis-eQTL (default: 1e6 bp)
- `pvalue_threshold_cis`: Significance threshold for cis-eQTLs
- `pvalue_threshold_trans`: Significance threshold for trans-eQTLs
- `error_covariance`: Covariance structure for error terms
- `use_model`: Statistical model (default: modelLINEAR)

**Covariates Included**:
- Age at index
- Gender
- Genetic ancestry (PC1-PC3)
- Methylation (gene-level beta values)
- CNV (gene-level copy number)

---

## Utility Scripts

### load_libraries.R

**Purpose**: Centralized library loading with dependency management

**Features**:
- Checks for installed packages
- Installs missing packages automatically
- Loads all required libraries
- Handles both CRAN and Bioconductor packages
- Provides version information

**Usage**:
```r
source('utils/load_libraries.R')
```

---

### data_utilities.R

**Purpose**: Common data manipulation and processing functions

**Key Functions**:
- `intersect_samples()` - Finds common samples across datasets
- `standardize_sample_ids()` - Normalizes TCGA sample IDs
- `impute_covariates()` - Handles missing values in covariates
- `create_covariate_matrix()` - Assembles covariate matrix from multiple sources

---

### plot_utilities.R

**Purpose**: Visualization helper functions

**Key Functions**:
- `plot_pca()` - Creates PCA scatter plots
- `plot_admixture()` - Plots admixture proportions
- `plot_manhattan()` - Creates Manhattan plots for eQTL results
- `plot_qq()` - Q-Q plots for p-value distributions
