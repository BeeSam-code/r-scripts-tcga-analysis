# TCGA R Scripts - Dependencies

## Core Data Manipulation & Processing

```r
library(dplyr)          # >= 1.0.0  - Data manipulation
library(readr)          # >= 2.0.0  - Fast file reading
library(tidyverse)      # >= 1.3.0  - Collection of data science packages
library(tidyr)          # >= 1.1.0  - Data tidying
library(data.table)     # >= 1.13.0 - Fast data frame alternative
library(readxl)         # >= 1.3.0  - Read Excel files
```

## Bioconductor Packages

```r
library(BiocManager)    # Bioconductor package manager
library(TCGAutils)      # >= 1.6.0  - TCGA data utilities
library(Biobase)        # >= 2.48.0 - Basic classes and methods
library(DESeq2)         # >= 1.28.0 - Differential expression analysis
library(edgeR)          # >= 3.30.0 - RNA-seq analysis
library(EnsDb.Hsapiens.v86)         # Ensembl annotation database
library(AnnotationDbi)  # >= 1.50.0 - Annotation database interface
library(GenomicFeatures) # >= 1.40.0 - Genomic feature extraction
```

## SNP & Genotype Analysis

```r
library(SNPRelate)      # >= 1.20.0 - SNP genome-wide association studies
library(genio)          # >= 1.1.0  - Genetics I/O utilities
library(SeqArray)       # >= 1.28.0 - Sequence array format
library(gdsfmt)         # >= 1.24.0 - GDS file format
library(trio)           # >= 3.20.0 - Family-based association analysis
```

## eQTL Analysis

```r
library(MatrixEQTL)     # >= 2.1.0  - Matrix eQTL for large-scale eQTL analysis
library(MASS)           # >= 7.3.0  - Modern applied statistics
```

## Normalization & Statistics

```r
library(preprocessCore) # >= 1.50.0 - Preprocessing and normalization
library(RNOmni)         # >= 1.0.0  - Rank-based inverse normal transformation
library(devtools)       # >= 2.3.0  - Development tools
```

## Visualization

```r
library(ggplot2)        # >= 3.3.0  - Grammar of graphics
library(magrittr)       # >= 2.0.0  - Pipe operators
```

## Utilities

```r
library(fs)             # >= 1.4.0  - File system operations
library(stringr)        # >= 1.4.0  - String manipulation
library(gdata)          # >= 2.18.0 - Data manipulation utilities
```

## Installation

### CRAN packages
```r
install.packages(c(
  'dplyr', 'readr', 'tidyverse', 'tidyr', 'data.table', 'readxl',
  'ggplot2', 'magrittr', 'fs', 'stringr', 'gdata', 'devtools', 'MASS'
))
```

### Bioconductor packages
```r
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install(c(
  'TCGAutils', 'Biobase', 'DESeq2', 'edgeR', 
  'EnsDb.Hsapiens.v86', 'AnnotationDbi', 'GenomicFeatures',
  'preprocessCore', 'RNOmni'
))
```

### SNPRelate and related packages
```r
BiocManager::install(c(
  'SNPRelate', 'SeqArray', 'gdsfmt', 'trio'
))
```

### Other packages
```r
install.packages('genio')
BiocManager::install('MatrixEQTL')
```

## System Requirements

- R >= 4.0.0
- Free disk space: ~5GB for typical TCGA datasets
- RAM: 16GB+ recommended for large-scale eQTL analysis

## Quick Installation Script

Use the `utils/load_libraries.R` script which handles installation and loading of all dependencies:

```r
source('utils/load_libraries.R')
```
