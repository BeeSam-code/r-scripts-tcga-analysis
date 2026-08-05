# TCGA Analysis R Scripts Repository

A clean, well-organized repository of R scripts for TCGA data analysis including gene expression, copy number variation (CNV), methylation, SNP/ancestry analysis, and eQTL identification.

## Repository Structure

```
r-scripts-tcga-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── CONTRIBUTING.md
├── data-processing/
│   ├── 01_rnaseq_data_processing.R
│   ├── 02_cnv_data_processing.R
│   └── 03_methylation_data_processing.R
├── analysis/
│   ├── 01_snp_ancestry_analysis.R
│   └── 02_eqtl_identification.R
├── utils/
│   ├── load_libraries.R
│   ├── data_utilities.R
│   └── plot_utilities.R
└── docs/
    ├── ANALYSIS_GUIDE.md
    ├── DEPENDENCIES.md
    └── SCRIPT_DESCRIPTIONS.md
```

## Scripts Overview

### Data Processing Scripts
- **01_rnaseq_data_processing.R**: Processes RNA-seq data from TCGA with normalization and filtering
- **02_cnv_data_processing.R**: Processes gene-level copy number variation data
- **03_methylation_data_processing.R**: Processes methylation array data with gene-level aggregation

### Analysis Scripts
- **01_snp_ancestry_analysis.R**: SNP analysis with ancestry inference using SNPRelate
- **02_eqtl_identification.R**: eQTL analysis using MatrixEQTL

## Quick Start

1. Clone: `git clone https://github.com/BeeSam-code/r-scripts-tcga-analysis.git`
2. Load dependencies: `source('utils/load_libraries.R')`
3. See [ANALYSIS_GUIDE.md](docs/ANALYSIS_GUIDE.md) for workflows

## Features

✅ Consistent naming conventions (snake_case)  
✅ Standard documentation headers with metadata  
✅ Modular functions extracted from procedural code  
✅ Removed duplicates and consolidated  
✅ Clean, consistent formatting  
✅ No hardcoded paths - parameterized inputs  
✅ Version-controlled and ready to use  

## License

MIT License - See LICENSE file
