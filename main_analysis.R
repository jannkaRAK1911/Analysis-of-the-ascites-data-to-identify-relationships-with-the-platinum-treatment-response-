install.packages("glmnet")
install.packages("randomForest")
install.packages("fastshap")
install.packages("iml")
install.packages("shapviz")
install.packages("ranger")
install.packages("xgboost")
install.packages("Matrix")
install.packages("SHAPforxgboost")
install.packages("caret")
install.packages("janitor")
library(readxl)
library(janitor)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(compositions)
library(glmnet)
library(randomForest)
library(fastshap)
library(iml)
library(shapviz)
library(ranger)
library(xgboost)
library(Matrix)
library(SHAPforxgboost)
library(glmnet)
library(caret)

################################################
######### Data loading + PREPROCESSING #########
################################################

dir.create(here('outputs', 'main_analysis'))

main_df <- read_excel(here("data", "main_dataset.xlsx"))
main_df <- main_df %>%
  filter(`Flow Cytometry` != 0)

# clr transformation of FC data
cell_colls <- c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
                "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells")
cells <- main_df[, cell_colls]
cells[cells == 0] <- 1e-6
clr_cells <- clr(cells)
main_df[, cell_colls] <- clr_cells

# log2 transformation of proteomic data 
prot_colls <- c("O00468",
                "O60240",                        
                "O75781", 
                "O95833",                        
                "P16050",                        
                "P17302",                       
                "P20292",                        
                "P23634",                        
                "P27487",                       
                "Q96CW1",                       
                "Q96Q06",                        
                "Q9BZV2",                       
                "Q9Y6N5")

main_df[, prot_colls] <- log2(main_df[, prot_colls] + 1)

# data scaling 
clinical_num_colls <- c("Ascites_Volume","Age_at_diagnosis")

num_colls <- c(clinical_num_colls, prot_colls, cell_colls)
main_df_scaled <- main_df
main_df_scaled[, num_colls] <- scale(main_df_scaled[, num_colls])

# Multinomial LASSO regression
cols_to_factor <- c("Diagnosis", "Stage_Clinical","Another_malignancy","NACHT", 
                    "PARPi", "BRCA_germinal", "BRCA_somatic", "P-response_primary_treatment",
                    "Flow Cytometry", "Organoids")
main_df_scaled_lasso <- main_df_scaled[, c(num_colls, cols_to_factor)]
main_df_scaled_lasso[is.na(main_df_scaled_lasso)] <- 0

# categorical cols -> factors
main_df_scaled_lasso <- main_df_scaled_lasso %>%
  mutate(across(all_of(cols_to_factor), as.factor))

main_df_scaled_lasso <- main_df_scaled_lasso %>%
  clean_names()
#y_lasso <- main_df_scaled_lasso$`P-response_primary_treatment`

bad_cols <- names(which(sapply(main_df_scaled_lasso, function(x)
  is.factor(x) && nlevels(x) < 2
)))

bad_cols

main_df_scaled_lasso <- main_df_scaled_lasso[, !names(main_df_scaled_lasso) %in% bad_cols]

#x_lasso <- model.matrix(`P-response_primary_treatment` ~ ., data = main_df_scaled_lasso)[, -1]

x <- main_df_scaled_lasso %>%
  select(-p_response_primary_treatment) %>%
  as.matrix()

y <- main_df_scaled_lasso$p_response_primary_treatment

cv_model <- cv.glmnet(
  x, y,
  family = "multinomial",  # 👈 kľúčové
  alpha = 1
)

# LASSO results
coef_lasso <- coef(cv_model, s = "lambda.min")

important_features <- unique(unlist(
  lapply(coef_lasso, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

important_features <- setdiff(
  important_features,
  "(Intercept)"
)

important_features

df_selected <- main_df_scaled_lasso %>%
  select(
    p_response_primary_treatment,
    all_of(important_features)
  )

## interakcie
safe_features <- paste0(
  "`",
  important_features,
  "`"
)

formula_inter <- as.formula(
  paste(
    "p_response_primary_treatment ~ (",
    paste(safe_features, collapse = " + "),
    ")^2"
  )
)

x_inter <- model.matrix(
  formula_inter,
  data = df_selected
)[,-1]

set.seed(123)

cv_model_inter <- cv.glmnet(
  x_inter,
  y,
  family = "multinomial",
  alpha = 1
)

coef_inter <- coef(
  cv_model_inter,
  s = "lambda.min"
)

selected_terms <- unique(unlist(
  lapply(coef_inter, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

selected_terms

interaction_terms <- selected_terms[
  grepl(":", selected_terms)
]

interaction_terms

coef_inter


## PCA
pca_data <- main_df_scaled_lasso %>%
  select(-p_response_primary_treatment)

pca_res <- prcomp(
  pca_data,
  center = TRUE,
  scale. = FALSE
)






























############################# New approach #################################


main_df_scaled_lasso$p_response_primary_treatment <- as.factor(main_df_scaled_lasso$p_response_primary_treatment
                                                              )

levels(main_df_scaled_lasso$p_response_primary_treatment)

x <- main_df_scaled_lasso %>%
  select(-p_response_primary_treatment) %>%
  as.matrix()

y <- main_df_scaled_lasso$p_response_primary_treatment

set.seed(123)

cv_fit <- cv.glmnet(
  x,
  y,
  family = "multinomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "class"
)

coef_list <- coef(cv_fit, s = "lambda.min")
names(coef_list)


important_features <- unique(unlist(
  lapply(coef_list, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

important_features <- setdiff(
  important_features,
  "(Intercept)"
)

important_features



df_selected <- main_df_scaled_lasso %>%
  select(all_of(c("p_response_primary_treatment", important_features)))


formula_inter <- as.formula(
  paste(
    "p_response_primary_treatment ~ (",
    paste(important_features, collapse = " + "),
    ")^2"
  )
)

formula_inter


x_inter <- model.matrix(
  formula_inter,
  data = df_selected
)[,-1]

y_inter <- df_selected$p_response_primary_treatment


set.seed(123)

cv_inter <- cv.glmnet(
  x_inter,
  y_inter,
  family = "multinomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "class"
)

coef_inter <- coef(
  cv_inter,
  s = "lambda.min"
)

selected_terms <- unique(unlist(
  lapply(coef_inter, function(x) {
    rownames(x)[x[,1] != 0]
  })
))

selected_terms


coef_inter$`1`
coef_inter$`3`
coef_list$`3`
