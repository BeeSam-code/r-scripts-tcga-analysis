#' ============================================================================
#' Utility: load_libraries.R
#' Purpose: Centralized library loading with automatic dependency management
#' Author: Sam Boudeau
#' Date: 2024-01-15
#'
#' Description:
#'   This script loads all required packages for TCGA analysis. It checks
#'   for installed packages and installs missing dependencies automatically.
#'   Handles both CRAN and Bioconductor packages.
#'
#' Dependencies:
#'   - BiocManager (for Bioconductor packages)
#'
#' Usage:
#'   source('utils/load_libraries.R')
#'
#' ============================================================================

# Function to safely load libraries
load_library <- function(package_name, source = 'CRAN') {
  if (!require(package_name, character.only = TRUE)) {
    message(sprintf("Installing %s from %s...", package_name, source))
    
    if (source == 'Bioconductor') {
      if (!require('BiocManager', character.only = TRUE)) {
        install.packages('BiocManager')
      }
      BiocManager::install(package_name)
    } else {
      install.packages(package_name)
    }
    
    require(package_name, character.only = TRUE)
    message(sprintf("%s loaded successfully", package_name))
  }
}

# Load CRAN packages
CRAN_PACKAGES <- c(
  'dplyr', 'readr', 'tidyverse', 'tidyr', 'data.table', 'readxl',
  'ggplot2', 'magrittr', 'fs', 'stringr', 'gdata', 'devtools',
  'MASS'
)

message("Loading CRAN packages...")
for (pkg in CRAN_PACKAGES) {
  load_library(pkg, source = 'CRAN')
}

# Load Bioconductor packages
BIOC_PACKAGES <- c(
  'TCGAutils', 'Biobase', 'DESeq2', 'edgeR',
  'EnsDb.Hsapiens.v86', 'AnnotationDbi', 'GenomicFeatures',
  'preprocessCore', 'RNOmni', 'SNPRelate', 'SeqArray',
  'gdsfmt', 'trio', 'MatrixEQTL'
)

message("Loading Bioconductor packages...")
for (pkg in BIOC_PACKAGES) {
  load_library(pkg, source = 'Bioconductor')
}

# Load other packages
load_library('genio', source = 'CRAN')

message("\n=" %+% rep('=', 78))
message("All libraries loaded successfully!")
message("Ready for TCGA analysis...")
message("=" %+% rep('=', 78))
