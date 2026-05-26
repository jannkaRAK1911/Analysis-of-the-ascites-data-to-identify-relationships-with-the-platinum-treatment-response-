################################################################################
########################### Factors identification #############################
################################################################################

#################################################################
# 1. LIBRARIES
#################################################################

# Library installation if needed:
install.packages("janitor")
install.packages("glmnet")


library(readxl)       # Read data from Excel files (.xlsx)
library(janitor)      # Clean data column names (clean_names)
library(here)         # Secure and robust relative file path management
library(dplyr)        # Data manipulation (filter, select, mutate, ...)
library(compositions) # CLR (Centered Log-Ratio) transformation for compositional data
library(glmnet)       # Regularized regression (LASSO, Elastic Net, Ridge)

#################################################################
# 2. DATA LOADING & PREPROCESSING
#################################################################

# Create an output directory if it doesn't exist yet
output_dir <- here('outputs', 'main_analysis')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the main dataset from Excel
main_df <- read_excel(here("data", "main_dataset.xlsx"))

# Filter: Keep rows where Flow Cytometry values are not zero
main_df <- main_df %>%
  filter(`Flow Cytometry` != 0)

# --- CLR Transformation for Flow Cytometry (FC) Data ---
# Define columns representing cell types (compositional data - parts of a whole)
cell_colls <- c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
                "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells")

cells <- main_df[, cell_colls]

# Handle zero values: Replace 0 with a small pseudocount 
# because log transformation in CLR cannot handle exact zeros
cells[cells == 0] <- 1e-6

# Apply Centered Log-Ratio (CLR) transformation and assign back to the dataset
clr_cells <- clr(cells)
main_df[, cell_colls] <- clr_cells

# --- Log2 Transformation for Proteomic Data ---
# Define columns containing proteomic markers (Uniprot IDs)
prot_colls <- c("O00468", "O60240", "O75781", "O95833", "P16050", 
                "P17302", "P20292", "P23634", "P27487", "Q96CW1", 
                "Q96Q06", "Q9BZV2", "Q9Y6N5")

# Log2 transformation (adding +1 to avoid log(0) for undetected proteins)
main_df[, prot_colls] <- log2(main_df[, prot_colls] + 1)

# --- Data Scaling (Standardization) ---
# Define numerical clinical variables
clinical_num_colls <- c("Ascites_Volume", "Age_at_diagnosis")

# Combine all numerical columns that need scaling
num_colls <- c(clinical_num_colls, prot_colls, cell_colls)
main_df_scaled <- main_df

# Z-score standardization (mean = 0, standard deviation = 1)
main_df_scaled[, num_colls] <- scale(main_df_scaled[, num_colls])

#################################################################
# 3. DATA PREPARATION FOR MULTINOMIAL LASSO
#################################################################

# Define categorical variables to be converted into factors
cols_to_factor <- c("Diagnosis", "Stage_Clinical", "Another_malignancy", "NACHT", 
                    "BRCA_germinal", "BRCA_somatic", "P-response_primary_treatment")

# Select relevant columns for the model (numerical features + categorical targets/features)
main_df_scaled_lasso <- main_df_scaled[, c(num_colls, cols_to_factor)]

# Convert selected columns into factors
main_df_scaled_lasso <- main_df_scaled_lasso %>%
  mutate(across(all_of(cols_to_factor), as.factor))

# Clean column names (removes special characters, converts to snake_case)
main_df_scaled_lasso <- main_df_scaled_lasso %>%
  clean_names()

# Detect "bad" columns: Identify categorical variables with less than 2 levels
# (variables with only 1 constant value provide no information for the model)
bad_cols <- names(which(sapply(main_df_scaled_lasso, function(x)
  is.factor(x) && nlevels(x) < 2
)))

# Print removed variables (if any)
bad_cols

# Exclude single-level factors from the dataset
main_df_scaled_lasso <- main_df_scaled_lasso[, !names(main_df_scaled_lasso) %in% bad_cols]

#################################################################
# 4. BASELINE LASSO MODEL (FEATURE SELECTION)
#################################################################

# Create the predictor matrix (X) by dropping the target variable and converting to matrix
x <- main_df_scaled_lasso %>%
  select(-p_response_primary_treatment) %>%
  as.matrix()

# Define the target variable (Y)
y <- main_df_scaled_lasso$p_response_primary_treatment

# Run cross-validated LASSO regression (alpha = 1)
cv_model <- cv.glmnet(
  x, y,
  family = "multinomial",  # Used for multi-class classification (target has > 2 levels)
  alpha = 1
)

# Extract non-zero coefficients at the optimal lambda value (lambda.min)
coef_lasso <- coef(cv_model, s = "lambda.min")

# Identify important features kept by LASSO across all outcome classes
important_features <- unique(unlist(
  lapply(coef_lasso, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

# Remove the Intercept term from the list of selected features
important_features <- setdiff(important_features, "(Intercept)")

# Print selected key features
important_features

# --- EXPORT 1 & 2: Baseline Model Outputs ---
# Save the list of baseline important features as a text file
writeLines(important_features, file.path(output_dir, "important_features_base.txt"))
# Save the complete baseline coefficient object (as an R RDS file)
saveRDS(coef_lasso, file.path(output_dir, "coef_base_multinomial.rds"))


#################################################################
# 5. LASSO MODEL WITH INTERACTIONS
#################################################################

# Subset the dataset to keep only the target variable and selected features
df_selected <- main_df_scaled_lasso %>%
  select(p_response_primary_treatment, all_of(important_features))

# Wrap feature names in backticks for safety (in case they contain special characters)
safe_features <- paste0("`", important_features, "`")

# Construct the formula for second-order interactions (all A * B combinations)
formula_inter <- as.formula(
  paste(
    "p_response_primary_treatment ~ (",
    paste(safe_features, collapse = " + "),
    ")^2"
  )
)

# Generate the model matrix including interaction terms (dropping the first column, which is the intercept)
x_inter <- model.matrix(formula_inter, data = df_selected)[,-1]

# Set seed for cross-validation reproducibility
set.seed(123)

# Run cross-validated LASSO regression on features including interactions
cv_model_inter <- cv.glmnet(
  x_inter, y,
  family = "multinomial",
  alpha = 1
)

# Extract coefficients for the interaction model at the optimal lambda
coef_inter <- coef(cv_model_inter, s = "lambda.min")

# Collect all terms (main effects and interactions) kept by the model
selected_terms <- unique(unlist(
  lapply(coef_inter, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

# Print all selected terms
selected_terms

# Filter: Extract only the terms that represent an actual interaction (containing a colon `:`)
interaction_terms <- selected_terms[grepl(":", selected_terms)]

# Final output of chosen interaction terms and complete model coefficients
interaction_terms
coef_inter

# --- EXPORT 3 & 4: Interaction Model Outputs ---
# Save the list of active interaction terms as a text file
writeLines(interaction_terms, file.path(output_dir, "interaction_terms.txt"))
# Save the complete interaction model coefficient object (as an R RDS file)
saveRDS(coef_inter, file.path(output_dir, "coef_interaction_multinomial.rds"))

