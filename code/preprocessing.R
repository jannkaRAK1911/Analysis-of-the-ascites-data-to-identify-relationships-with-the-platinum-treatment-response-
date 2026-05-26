################################################################################
################ INITIAL PROCESSING – PROTEOMICS ###############################
################################################################################

################################################################################
# LIBRARIES 
################################################################################

# Load required libraries for:
# - gene symbol validation and correction
# - data manipulation
# - proteomics preprocessing
# - reading Excel files
library(HGNChelper)
library(here)
library(dplyr)
library(tidyr)
library(DEP)
library(readxl)

################################################################################
# DATA LOADING 
################################################################################

# Create output directory for processed data and QC files
dir.create(here("outputs", "data-processing"),
           recursive = TRUE,
           showWarnings = FALSE)

# Define input file and relevant sheet/metric
xlsx_path <- here("data/6160_PG master report_250106.xlsx")
sheet <- "Protein Groups"
metric <- "PG MaxLFQ Intensities - log2 Transformed"

# Load Excel sheet without predefined column names
df <- read_excel(xlsx_path,
                 sheet = sheet,
                 col_names = FALSE)

# Extract proper column names from the second row
colnames(df) <- as.character(df[2, ])

# Remove metadata/header rows
df <- df[-c(1, 2), ]

################################################################################
# SAMPLE SELECTION 
################################################################################

# Select relevant columns:
# - protein identifiers
# - intensity measurements
# - response annotation
#
# Different intensity types can be selected depending on analysis:
# - normalized MaxLFQ
# - non-normalized MaxLFQ
# - raw intensities

# MaxLFQ not normalized
df_55 <- df[, c(1:3, 370:424, 540)]

# MaxLFQ normalized
# df_55 <- df[, c(1:3, 426:480, 540)]

# Raw intensities
# df_55 <- df[, c(1:3, 255:310, 540)]

# Extract correct sample names from the third row
colnames(df_55) <- as.character(df_55[3, ])

# Remove remaining technical/header rows
df_55_final <- df_55[-c(1, 2, 3, 4), ]

# Rename gene name column for compatibility with DEP workflow
proteinGroups <- df_55_final %>%
  rename(Gene.Name = `Gene Name`)

# Save unfiltered protein table for QC/reference
write.csv(
  proteinGroups,
  here("outputs", "data-processing", "proteinGroups_cRAP.csv"),
  row.names = FALSE
)

################################################################################
# CONTAMINANT REMOVAL 
################################################################################

# Identify common contaminant proteins originating from cRAP database
# (MaxQuant labels these proteins with "cRAP" in accession names)
cRAP.proteins <- proteinGroups %>%
  filter(grepl("cRAP", Accession))

# Save contaminant list for quality control purposes
write.csv(
  cRAP.proteins,
  here("outputs", "data-processing", "proteinGroups_cRAP.csv"),
  row.names = FALSE
)

# Remove contaminant proteins from the dataset
d1 <- proteinGroups %>%
  filter(!grepl("cRAP", Accession))

################################################################################
# OPTIONAL FILTERING
################################################################################

# Remove keratin proteins if protein descriptions are available
# Keratins are common contaminants in proteomics experiments
if ("Description" %in% colnames(d1)) {

  d1 <- d1 %>%
    filter(!grepl("keratin",
                  Description,
                  ignore.case = TRUE))
}

################################################################################
# UNIQUE IDENTIFIER GENERATION
################################################################################

# Generate unique protein/gene identifiers required by DEP package
# This prevents duplicated gene symbol issues during downstream analysis
d1 <- make_unique(
  d1,
  "Accession",
  "Gene.Name",
  delim = ";"
)

################################################################################
# GENE SYMBOL HARMONIZATION 
################################################################################

# Load HGNC-approved gene symbol mapping table
load(here("data", "genenames_update_20220816.RData"))

# Validate and update gene symbols
update_d1 <- checkGeneSymbols(
  d1$Gene.Name,
  species = "human",
  map = genenames_newest,
  unmapped.as.na = FALSE
)

################################################################################
# MANUAL AMBIGUITY RESOLUTION 
################################################################################

# Create lookup table:
# old gene symbol -> approved HGNC symbol
gene_map <- setNames(genenames_newest$Approved.Symbol,
                     genenames_newest$Symbol)

# Resolve ambiguous gene annotations:
# - split multi-gene entries
# - map outdated symbols
# - retain the first approved symbol
resolve_gene <- function(x, gene_map) {

  parts <- unlist(strsplit(x,
                           "\\s*(///|;|\\?)\\s*"))

  parts_mapped <- sapply(parts, function(g) {

    if(g %in% names(gene_map)) {
      return(gene_map[g])
    } else {
      return(g)
    }
  })

  return(parts_mapped[1])
}

# Apply ambiguity resolution to all suggested symbols
update_d1$Suggested.Symbol <- sapply(
  update_d1$Suggested.Symbol,
  resolve_gene,
  gene_map
)

# Manual corrections for specific unresolved conflicts
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "EPRS1 /// QARS1"] <- "QARS1"

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "C12orf75 /// FSTL1"] <- "C12orf75"

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "LNPK /// NUSAP1"] <- "LNPK"

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "SARS1 /// SARS2"] <- "SARS1"

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol ==
                             "MPHOSPH6 /// PALS2"] <- "PALS2"

# Inspect duplicated suggested gene symbols
update_d1$Suggested.Symbol[
  duplicated(update_d1$Suggested.Symbol)
]

# Merge updated gene annotations back into the main dataset
d1 <- left_join(
  d1,
  update_d1,
  by = c("Gene.Name" = "x")
)

# Additional duplicate inspection after merging
d1$Suggested.Symbol[
  duplicated(d1$Suggested.Symbol)
]

duplicates <- d1 %>%
  filter(Gene.Name %in%
           d1$Gene.Name[duplicated(d1$Gene.Name)])

dup_genes <- d1 %>%
  group_by(Gene.Name) %>%
  filter(n() > 1) %>%
  arrange(Gene.Name)

# Remove fully duplicated rows
d1 <- d1[!duplicated(d1), ]

################################################################################
# FINAL CLEANUP
################################################################################

# Rename mislabeled sample
colnames(d1)[colnames(d1) == "245_NA"] <- "245_SE"

################################################################################
# FINAL EXPORT 
################################################################################

# Save cleaned and harmonized proteomics dataset
write.csv(
  d1,
  here("outputs", "data-processing", "updated_proteinGroups_SEC_only.csv"),
  row.names = FALSE
)

################################################################################
############################ QC SUMMARY ########################################
################################################################################

message("Initial proteins: ", nrow(proteinGroups))
message("Removed cRAP proteins: ", nrow(cRAP.proteins))
message("Final proteins after cleaning: ", nrow(d1))
