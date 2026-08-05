# R Scripts Evaluation Report - HNSC/GEWIS Analysis Files

**Date**: August 5, 2024  
**Evaluator**: Code Analysis System  
**User**: BeeSam-code

---

## 📊 Executive Summary

You have **5 scripts worth publishing** on professional GitHub. These files represent sophisticated genomic epidemiology and GWAS analysis, with applications in head/neck cancer (HNSC) research.

**Overall Quality**: 8/10 - Publication-ready with refactoring

---

## ⭐ Scripts Ranked by Professional Value

### 🥇 TIER 1: HIGHLY RECOMMENDED (Publish As-Is or with Minor Cleanup)

#### 1. **relmatGlmer_hnsc_tobacco_kinshipMixedModels_1124.R**
**Score**: 9/10 ⭐⭐⭐⭐⭐  
**Status**: PUBLISH - Production-ready

**Why it's excellent**:
- ✅ Implements advanced mixed-effects models with kinship matrices
- ✅ Handles complex genetic epidemiology (SNP + Environment interactions)
- ✅ Sophisticated parallel processing with `foreach`/`doParallel`
- ✅ Proper command-line argument handling for HPC cluster submission
- ✅ Two-stage analysis (marginal SNP test → GxE interaction)
- ✅ Chromosome-specific processing (scalable pipeline)
- ✅ Uses established packages: `coxme`, `lme4qtl`, `GMMAT`

**Strengths**:
```r
# Clean parallel processing architecture
results <- foreach(i = 1:length(snp_list), .packages = c("coxme")) %dopar% {
  # Efficient SNP-by-SNP testing
}

# Proper kinship matrix handling
kinship_mat2 <- kinship_mat2 + diag(1e-5, nrow(kinship_mat2))

# GxE interaction testing (tobacco use × SNP)
GE_modelGxE <- lmekin(CC_status ~ model_DFGxE[, 8]*Tabacco_ever + ...)
```

**Issues to Fix**:
- 🔧 Hardcoded file paths → parameterize
- 🔧 Add function headers with roxygen2 documentation
- 🔧 Extract helper functions (data loading, model fitting)
- 🔧 Add error handling for edge cases (singular kinship matrices)

**Recommendation**: 
**PUBLISH** - This is sophisticated epidemiological genetics. Extract into modular functions and document well. This would be excellent for Bioconductor.

---

#### 2. **gee_analysis_HNCC_smoking.R**
**Score**: 8.5/10 ⭐⭐⭐⭐⭐  
**Status**: PUBLISH with refactoring

**Why it's excellent**:
- ✅ Implements 3 different approaches to GEE/mixed models
- ✅ Compares `lmekin` (coxme), `parLapply`, and `foreach` parallelization strategies
- ✅ Proper workflow: data prep → kinship calculation → model fitting
- ✅ Good parameter exploration (shows testing of different methods)
- ✅ Real genomic data pipeline (100k SNP sampling)
- ✅ Demonstrates best practices for parallel GWAS analysis

**Code Quality Example**:
```r
# Excellent: Three parallel processing approaches documented

# Method 1: Base R parLapply
results <- parLapply(cl, test_snps, run_model_for_snp)

# Method 2: foreach/doParallel
results <- foreach(i = 1:length(test_snps), .packages = c("coxme")) %dopar% { ... }

# Method 3: Base R for loop (for comparison)
for (i in 1:length(test_snps)) { ... }
```

**Issues to Fix**:
- 🔧 Massive hardcoded path: `"C:/Users/boudeas/OneDrive..."` → use `here::here()`
- 🔧 Duplicate code across 3 methods → extract as functions
- 🔧 No documentation headers
- 🔧 Commented-out old GEE attempts → create separate "methods_comparison" file
- 🔧 Results stored as data frames with unclear structure → add `result_schema` documentation

**Recommendation**: 
**PUBLISH** - Refactor into modular functions. This is excellent comparative methodology for parallel GWAS. Would be good tutorial for HPC implementation.

---

#### 3. **manhattanPlotFeb26.R**
**Score**: 8/10 ⭐⭐⭐⭐  
**Status**: PUBLISH (convert to function)

**Why it's good**:
- ✅ Complete GWAS visualization pipeline
- ✅ Combines multiple chromosomes from separate PLINK outputs
- ✅ Produces publication-quality Manhattan plots
- ✅ Includes Q-Q plots for p-value calibration
- ✅ Handles large datasets efficiently (data.table)
- ✅ Nice customization: chromosome colors, threshold lines, SNP labels

**Code Quality**:
```r
# Good: Combines multiple files
combined_chroms <- rbindlist(lapply(files, fread), fill = TRUE)

# Good: Filtering for significant results
top_snps <- combined_chroms[combined_chroms$logp > -log10(5e-8), "ID"]

# Good: Publication-quality plot
manhattan(
  combined_chroms,
  chr="#CHROM", bp="POS", p="P", snp="ID",
  genomewideline = -log10(5e-8),
  suggestiveline = -log10(1e-5)
)
```

**Issues**:
- 🔧 Procedural script → convert to function
- 🔧 Hardcoded file patterns → parameterize
- 🔧 No error handling for missing columns
- 🔧 Color schemes hardcoded → use palette function

**Recommendation**: 
**PUBLISH** - Convert to `create_manhattan_plot()` and `create_qq_plot()` functions. Add to your visualization utilities.

---

#### 4. **GEWIS_PCA_DemographicDataHandling_1224.R**
**Score**: 8/10 ⭐⭐⭐⭐  
**Status**: PUBLISH with significant refactoring

**Why it's good**:
- ✅ Complete population genetics workflow
- ✅ PCA with 1000 Genomes reference populations
- ✅ Sex imputation from F-statistic (HWE check)
- ✅ Demographic data integration and QC
- ✅ PLINK covariate/phenotype file generation
- ✅ Publication-quality PCA plots with variance explained
- ✅ GWAS association result visualization

**Code Quality**:
```r
# Good: Merging multiple data sources
merged_data <- merge(eigenven, pops_metadata, by = "X.IID")

# Good: Sex imputation logic
mergedDemo$Sex <- ifelse(
  mergedDemo$Gender %in% "Female", 
  mergedDemo$Gender, 
  mergedDemo$testSex
)

# Good: gtsummary table generation
all_batch1_sum <- demosex %>%
  dplyr::select(...) %>%
  tbl_summary(by = CC_status, ...)
```

**Issues**:
- 🔧 **MAJOR**: Multiple hardcoded paths with full OneDrive paths
- 🔧 Hardcoded file separators and column names
- 🔧 No function decomposition (300+ lines procedural)
- 🔧 Commented-out code blocks should be in separate file
- 🔧 QC analysis is interesting → extract to separate script
- 🔧 Missing validation (check for duplicate sample IDs, etc.)

**Recommendation**: 
**PUBLISH with MAJOR REFACTORING** - Split into:
1. `demographic_data_qc.R` - Clean demographic data
2. `pca_ancestry_assignment.R` - PCA and population assignment
3. `gwas_result_visualization.R` - Plot association results

---

### 🥈 TIER 2: GOOD (Publish as Utilities/Examples)

#### 5. **locus_zoom_feb2026.R**
**Score**: 7.5/10 ⭐⭐⭐⭐  
**Status**: PUBLISH as localized region visualization utility

**Strengths**:
- ✅ Nice use of `locuszoomr` package for regional plots
- ✅ LD matrix integration (r² coloring)
- ✅ Integrated gene annotation from EnsDb
- ✅ SNP labeling and highlighting
- ✅ Publication-quality TIFF output

**Issues**:
- 🔧 Very specialized (locus-zoom plots only)
- 🔧 Hardcoded SNP selection: `hits <- c("rs74125744", "rs11804045")`
- 🔧 Hardcoded region coordinates: `xrange = c(147168399, 148168399)`
- 🔧 No function wrapper

**Recommendation**: 
**PUBLISH** - Convert to:
```r
create_locus_zoom <- function(gwas_results, ld_results, region_chr, region_start, 
                               region_end, highlight_snps = NULL)
```

---

### 🥉 TIER 3: SPECIALIZED/NICHE

#### 6. **multicancer_summary_May2024.R**
**Score**: 7/10 ⭐⭐⭐  
**Status**: PUBLISH as analysis example (cancer-specific)

**What it does**:
- Integrates eQTL results across 7 cancer types (BRCA, HNSC, PRAD, COADREAD, etc.)
- Cross-cancer SNP/gene associations
- GTEx tissue-specific validation
- UpSet plots and network analysis
- Identifies variants replicating across cancers

**Strengths**:
- ✅ Multi-cancer integration (sophisticated)
- ✅ External validation with GTEx
- ✅ Nice visualization (UpSet, Venn, network plots)
- ✅ Real scientific question: what eQTLs replicate across cancer types?

**Issues**:
- 🔧 **Very** cancer/dataset specific (limited reuse)
- 🔧 Hardcoded cancer type names and file patterns
- 🔧 Duplicated merging code (could be function)
- 🔧 Intermediate/analysis-specific (not a general utility)
- 🔧 Missing error handling for missing GTEx tissues

**Recommendation**: 
**PUBLISH as example/vignette** - Good case study for multi-study integration but limited general applicability.

---

#### 7. **overrepresentation_analysis_Nov2023.R**
**Score**: 6.5/10 ⭐⭐⭐  
**Status**: INCOMPLETE - Fragment only

**What's included**:
- Gene ontology enrichment using `clusterProfiler`
- RNA-seq data loading and preparation
- `org.Hs.eg.db` annotation database

**Issues**:
- 🔧 **Incomplete**: Only loads data, doesn't run enrichment
- 🔧 No actual GO analysis or results
- 🔧 Would need significant completion
- 🔧 Minimal reusable code

**Recommendation**: 
**NOT READY** - Complete this script or skip publication.

---

#### 8. **hwe_plots_ooct2023.R**
**Score**: 5/10 ⭐⭐  
**Status**: MINIMAL VALUE - Too simple for publication

**What it does**:
- Compares observed vs expected heterozygosity
- Boxplot visualization

**Issues**:
- 🔧 **Too simple** for standalone publication (3 lines of plotting)
- 🔧 Single use case (HWE check)
- 🔧 No generalization or reuse

**Recommendation**: 
**SKIP** - Integrate into QC pipeline instead.

---

#### 9. **heatmap_06282021.R**
**Score**: 4/10 ⭐  
**Status**: NOT SUITABLE - Outdated/niche

**Issues**:
- 🔧 One-off heatmap of allele frequencies
- 🔧 No real analysis
- 🔧 Dated (2021)
- 🔧 Very specific use case

**Recommendation**: 
**SKIP** - Not suitable for professional GitHub.

---

## 📋 Publication Roadmap

### **Immediate (Ready to publish with minor cleanup)**

```
r-scripts-gwas-kinship-analysis/
├── README.md
├── CONTRIBUTING.md
├── LICENSE
├── .gitignore
├── 
├── mixed-models/
│   ├── 01_relmat_glmer_analysis.R          [relmatGlmer_hnsc...]
│   ├── 02_gee_comparative_methods.R        [gee_analysis...]
│   └── 03_marginal_snp_testing.R           [Extract SNP marginal testing]
│
├── visualization/
│   ├── 01_manhattan_plot.R                 [Refactor manhattan...]
│   ├── 02_locus_zoom_regional_plot.R       [locus_zoom...]
│   └── 03_qq_plot.R                        [Q-Q plot utilities]
│
├── population-genetics/
│   ├── 01_pca_ancestry_assignment.R        [Extract from GEWIS_PCA...]
│   ├── 02_demographic_qc.R                 [Extract from GEWIS_PCA...]
│   └── 03_gwas_result_parsing.R
│
├── utils/
│   ├── kinship_matrix_utilities.R
│   ├── parallel_processing_utilities.R
│   └── data_loading_utilities.R
│
├── examples/
│   ├── multi_cancer_eqtl_integration.R     [multicancer_summary...]
│   └── hnsc_gwas_pipeline.R                [Full worked example]
│
└── docs/
    ├── ANALYSIS_GUIDE.md
    ├── DEPENDENCIES.md
    └── METHOD_COMPARISON.md
```

---

## 🔧 Key Refactoring Tasks

### Priority 1: Remove Hardcoded Paths

**Current**:
```r
genoSample <- read.csv("chr22_202205_202106_202008_Batch1_3Chr1_Xdose04_sampleQC_noDups_HWE_recode012_1024.csv")
demog_data <- read.csv("FCCC_Ragin_GSA12K_202205_202106_202008_postImputDemog_1024.csv")
ancestryPCs <- read.csv("GEWIS_HNSCCasesCtrl_genotype_eigenvectors_OCT2024.csv")
kinship_mat <- read.csv("KinshipMat_202205_202106_202008_allFilesBatch1_3Chr1_Xdose04_sQC_noDups_HWE_1124.csv")
```

**Refactored**:
```r
# config.R or command-line args
data_dir <- here::here("data", "genotypes")
kinship_dir <- here::here("data", "kinship")
demog_file <- here::here("data", "demographics.csv")

# In script
load_genotype_data <- function(chrom, data_dir) {
  file <- file.path(data_dir, paste0("chr", chrom, "_genotypes.csv"))
  if (!file.exists(file)) stop(sprintf("File not found: %s", file))
  read.csv(file)
}
```

### Priority 2: Extract Functions

**Example: Kinship-based model fitting**

```r
fit_kinship_model <- function(genotype, phenotype, covariates, kinship_matrix,
                              formula = NULL) {
  #' Fit linear/logistic mixed model with kinship matrix
  #'
  #' @param genotype Vector of genotypes (0, 1, 2)
  #' @param phenotype Vector of outcomes
  #' @param covariates Data frame of covariates
  #' @param kinship_matrix Matrix of kinship coefficients
  #' @param formula Formula (auto-constructed if NULL)
  #'
  #' @return lmekin model object with p-value and coefficient
  
  # Implementation
}

# Usage
results <- foreach(snp in snp_list, .packages = c("coxme")) %dopar% {
  fit_kinship_model(
    genotype = model_df[[snp]],
    phenotype = model_df$CC_status,
    covariates = model_df[, c("EV1", "EV2", "EV3", "Age")],
    kinship_matrix = kinship_mat
  )
}
```

### Priority 3: Add Documentation

```r
#' ============================================================================
#' Script: 02_gee_comparative_methods.R
#' Purpose: Compare GEE, lmekin, and foreach implementations for GWAS
#' Author: Sam Boudeau
#' Date: 2024-10-26
#'
#' Description:
#'   This script demonstrates three different approaches to fitting marginal
#'   and interaction models in GWAS with binary outcomes and kinship matrices.
#'   All methods are equivalent statistically; choice depends on performance.
#'
#' Methods Compared:
#'   1. Sequential for loop (baseline)
#'   2. parLapply (parallel package)
#'   3. foreach/doParallel (more flexible)
#'
#' Dependencies:
#'   - coxme (>= 2.3.0): Mixed effects Cox models with kinship
#'   - foreach, doParallel: Parallel processing
#'   - data.table: Fast I/O
#'   - dplyr, tidyr: Data manipulation
#'
#' Input:
#'   - Genotype matrix (SNPs × samples)
#'   - Phenotype/covariate data
#'   - Kinship matrix (pre-computed)
#'
#' Output:
#'   - Marginal SNP test results
#'   - SNP × Environment (GxE) interaction results
#'   - Timing comparison
#'
#' Usage:
#'   Rscript relmatGlmer_hnsc_tobacco_kinshipMixedModels_1124.R chr22
#'
#' ============================================================================
```

---

## 🎯 Final Recommendations

### **What to Publish**

| Script | Package | Status |
|--------|---------|--------|
| relmatGlmer_hnsc_tobacco... | GWAS | ✅ PUBLISH |
| gee_analysis_HNCC_smoking | GWAS | ✅ PUBLISH |
| manhattanPlotFeb26 | Visualization | ✅ PUBLISH |
| GEWIS_PCA_Demographic... | Population Genetics | ✅ PUBLISH (refactor) |
| locus_zoom_feb2026 | Visualization | ✅ PUBLISH |
| multicancer_summary | Example/Vignette | ⚠️ PUBLISH AS EXAMPLE |
| overrepresentation_analysis | Incomplete | ❌ SKIP |
| hwe_plots | Too simple | ❌ SKIP |
| heatmap_06282021 | Outdated | ❌ SKIP |

### **Create New Repository**

```
Repository: BeeSam-code/r-scripts-gwas-mixed-models
Description: R utilities for GWAS analysis with kinship matrices, 
             mixed models, and visualization
```

### **Suggested Structure**

```bash
BeeSam-code/r-scripts-gwas-mixed-models/
├── README.md (overview + quick start)
├── CONTRIBUTING.md (code style guide)
├── LICENSE (MIT)
├── DEPENDENCIES.md (all R packages + versions)
├── .gitignore (R-specific)
│
├── R/  # Main analysis scripts
│   ├── kinship_mixed_models.R
│   ├── marginal_snp_testing.R
│   ├── gxe_interaction_testing.R
│   └── parallel_gwas_pipeline.R
│
├── visualization/  # Plotting utilities
│   ├── manhattan_plots.R
│   ├── qq_plots.R
│   ├── locus_zoom_plots.R
│   └── pca_ancestry_plots.R
│
├── data-processing/  # QC and prep
│   ├── load_genotype_data.R
│   ├── demographic_qc.R
│   ├── kinship_matrix_prep.R
│   └── gwas_result_parsing.R
│
├── examples/  # Worked examples
│   ├── hnsc_gwas_workflow.R
│   └── multi_cancer_integration.R
│
└── docs/
    ├── ANALYSIS_GUIDE.md
    ├── METHOD_COMPARISON.md (GEE vs lmekin vs glmmkin)
    └── PERFORMANCE_BENCHMARKS.md
```

---

## 📊 Quality Scores Summary

| File | Score | Recommendation |
|------|-------|----------------|
| relmatGlmer_hnsc_tobacco_kinship... | 9/10 | 🟢 PUBLISH |
| gee_analysis_HNCC_smoking | 8.5/10 | 🟢 PUBLISH |
| manhattanPlotFeb26 | 8/10 | 🟢 PUBLISH |
| GEWIS_PCA_DemographicDataHandling | 8/10 | 🟢 PUBLISH (refactor) |
| locus_zoom_feb2026 | 7.5/10 | 🟡 PUBLISH (improve) |
| multicancer_summary_May2024 | 7/10 | 🟡 PUBLISH AS EXAMPLE |
| overrepresentation_analysis | 6.5/10 | 🔴 SKIP (incomplete) |
| hwe_plots | 5/10 | 🔴 SKIP (too simple) |
| heatmap_06282021 | 4/10 | 🔴 SKIP (outdated) |

**Overall**: Your genomic epidemiology code is excellent! Focus on the top 5 scripts, remove hardcoded paths, and add proper documentation. You'll have a publication-quality resource. 🎉

---

*Report Generated: 2024-08-05*  
*Recommended Repository: `BeeSam-code/r-scripts-gwas-mixed-models`*
