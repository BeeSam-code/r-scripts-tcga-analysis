# Contributing Guidelines

## Code Style

### Naming Conventions
- **Functions**: `snake_case` (e.g., `load_rnaseq_data`, `calculate_normalization`)
- **Variables**: `snake_case` (e.g., `gene_counts`, `sample_metadata`)
- **Constants**: `SCREAMING_SNAKE_CASE` (e.g., `MIN_EXPRESSION_THRESHOLD`)
- **Files**: `snake_case.R` with numbered prefix for processing order

### Formatting
- **Indentation**: 2 spaces
- **Line length**: Max 80 characters
- **Comments**: Use `#` with descriptive text explaining WHY, not WHAT

## Script Headers

All scripts must include:

```r
#' ============================================================================
#' Script: script_name.R
#' Purpose: Brief description
#' Author: Your name
#' Date: YYYY-MM-DD
#'
#' Description:
#'   Detailed description of functionality
#'
#' Dependencies:
#'   - package1 (>= version)
#'   - package2 (>= version)
#'
#' Input:
#'   - file.csv: Description
#'
#' Output:
#'   - output.csv: Description
#' ============================================================================
```

## Before Submitting

1. Test your code
2. Follow naming conventions
3. Remove commented code
4. Add documentation headers
5. Check formatting consistency

## Submitting

1. Create a feature branch
2. Commit with clear messages
3. Submit a Pull Request
