# Analysis of the ascites data to identify relationships with the platinum treatment response
Analytical workflows in R for the identification of factors associated with specific platinum-based chemotherapy response in ovarian cancer. The scripts integrate clinical data, extracellular vesicle proteomics, and flow cytometry data obtained from ascites. Using multinomial LASSO regression, it identifies key factors and their interactions associated with treatment outcome.

## Pipeline Overview
The pipeline is divided into three main functional scripts:
1. **Data Preprocessing & Gene Harmonization** (`preprocessing.R`)
2. **Compositional Data Analysis & Multi-platform Validation** (`cell_composition_analysis.R`)
3. **Feature Selection & Interaction Modeling** (`factors_identification.R`)

---

## Repository Structure & Script Descriptions

### 1. `preprocessing.R` (Initial Processing – Proteomics)
This script handles the raw proteomics discovery data (MaxLFQ intensities) and performs rigorous quality control (QC) and biological identifier standardization.
* **Data Loading & Cleaning:** Imports specific metrics from multi-sheet Excel reports and extracts proper sample mappings.
* **Contaminant Filtering:** Screens and removes common technical contaminants (e.g., from the cRAP database, keratins).
* **Gene Symbol Harmonization:** Leverages the `HGNChelper` package to validate and update outdated or ambiguous gene symbols to the latest approved HGNC nomenclature. It includes custom manual ambiguity resolution for multi-gene clusters (e.g., *SEPTIN2/6*, *EPRS1/QARS1*).
* **Outputs:** Generates a filtered, deduplicated, and harmonized protein intensity table (`updated_proteinGroups_SEC_only.csv`).

### 2. `cell_composition_analysis.R` (Cell Composition Analysis)
A dedicated Compositional Data Analysis (CoDA) and validation script that correlates EV-derived signatures with physical cell-type benchmarks.
* **Marker Targeting:** Filters EV profiles using specific cell-type marker genes derived from literature (Izar et al.) to estimate relative tissue sub-compositions (Malignant, Fibroblast, Macrophage).
* **Multi-platform Cross-Validation:** Uses **Repeated Measures Correlation (`rmcorr`)** and Spearman rank evaluations to statistically prove how well shedding EV protein signatures reflect true underlying Flow Cytometry (FC) percentages.
* **Multivariate Statistics:** Runs **PERMANOVA (`adonis2`)** and **Principal Component Analysis (PCA)** on Aitchison geometry spaces to see if multivariate cellular and vesicular architectures shift based on therapeutic response groups (*RT*, *RP*, *SE*).
* **Visualizations:** Generates high-quality PDF publications containing stacked bar charts, stratified treatment group layouts (built with `patchwork`), and pairwise Wilcoxon boxplots on CLR-transformed ratios.

### 3. `factors_identification.R` (Multinomial LASSO Regression)
This script builds a predictive model to identify which molecular (proteomic) and cellular components are tightly associated with patient clinical characteristics and treatment responses.
* **Compositional & Log Transformations:** Applies Centered Log-Ratio (CLR) transformations to cellular data and $Log_2$ transformations to proteomic data to ensure statistical validity.
* **Baseline Feature Selection:** Runs a cross-validated (`cv.glmnet`) **Multinomial LASSO regression** to isolate zero-variance categorical factors and rank key predictors at the optimal penalty level (`lambda.min`).
* **Interaction Modeling:** Constructs a second-order interaction matrix ($(A \times B)^2$) using the selected baseline features to discover complex cross-talks/co-dependencies.
* **Outputs:** Exports raw lists of key drivers (`important_features_base.txt`, `interaction_terms.txt`) and complete RDS model coefficient objects directly into the `outputs/main_analysis/` directory.

---

## Installation & Required Packages

To run this pipeline, you will need **R (version 4.5.2)**. You can install all required dependencies by running the following snippet inside your R console:

```R
required_packages <- c(
  "readxl", "janitor", "here", "dplyr", "tidyr", 
  "compositions", "glmnet", "patchwork",
  "ggplot2", "rmcorr", "ggprism", 
  "vegan", "ggpubr", "stringr", "tibble"
)

install.packages(setdiff(required_packages, installed.packages()[,"Package"]))
```
Packages HGNChelper and DEP can be installed using the following commands:

```R
if (!require("BiocManager"))
    install.packages("BiocManager")

BiocManager::install(c("HGNChelper", "DEP"))
