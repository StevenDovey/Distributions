# Distributions
R pipeline for NZ Pinus radiata forest structure analysis. Standardizes regional inventory datasets, applies cohort-specific 0–1 normalization to isolate size shifts, and fits crash-proof competitive distributions (Weibull, Gamma, Lognormal) via parallel processing. Exports AIC summaries and parameter metrics for trend modelling.



```markdown
# NZ Pinus radiata Forest Structure Analysis Pipeline

An optimized R pipeline to analyze structural dynamics in New Zealand *Pinus radiata* plantations across multiple decades and regions. This repository standardizes raw inventory data, applies cohort-specific normalization, evaluates distribution models via AIC, and generates high-throughput visualizations using parallel processing.

## Features

- **Data Standardization:** Maps raw geographic logging codes to official NZ regional classifications.
- **Proportional Scaling (0 to 1):** Scales tree DBH and height relative to cohort maximums to isolate structural profile shifts from absolute size increases.
- **Crash-Proof MLE Modeling:** Uses `tryCatch` to safely fit distribution parameters across sparse or skewed early-age classes.
- **Competitive Distribution Fit:** Runs competitive MLE fits across Weibull, Gamma, and Lognormal distributions to export an AIC selection summary.
- **Parallel Plotting Engine:** Leverages `foreach` and `doParallel` to distribute plot generation across multiple CPU cores.

## Project Structure

```text
├── forestry_data.csv                    # Input inventory dataset (untracked)
├── run_pipeline.R                       # Main execution script
├── best_fit_distributions_summary.csv   # Exported AIC model selection table
└── weibull_parameters_trend_analysis.csv # Exported shape/scale metrics for ANOVA/Regression

```

## Getting Started

### Prerequisites

```R
install.packages(c("tidyverse", "fitdistrplus", "foreach", "doParallel"))

```

### Execution

```R
source("run_pipeline.R")

```

The script automatically detects available CPU cores, spins up a parallel cluster, writes the statistical summary tables, and saves the rendering cohorts into structured age folders (`Age_X/`).

## Outputs

1. **`best_fit_distributions_summary.csv`**: Used for Chi-Square tests of independence or AIC delta evaluations to track distribution shifts across forest management eras.
2. **`weibull_parameters_trend_analysis.csv`**: Tracks Weibull shape ($\beta$) and scale ($\eta$) parameters for Two-Way ANOVAs or Linear Mixed Models to test variations in forest uniformity over time.
