################################################################################
################ INITIAL PROCESSING – SEC-only PROTEOMICS ######################
################################################################################


library(HGNChelper)
library(here)
library(dplyr)
library(tidyr)
library(DEP)
library(readxl)
library(UpSetR)

dir.create(here("outputs", "data-processing"), recursive = TRUE, showWarnings = FALSE)

# --- INPUT ---
xlsx_path <- here("data/6160_PG master report_250106.xlsx")  
sheet <- "Protein Groups"
metric <- "PG MaxLFQ Intensities - log2 Transformed"

# Načítaj data
df <- read_excel(xlsx_path, sheet = sheet, col_names = FALSE)

colnames(df) <- as.character(df[2, ])

df <- df[-c(1, 2), ]

#maxLFQ normalized
#df_55 <- df[,c(1:3, 426:480, 540)]
#maxLFQ raw
df_55 <- df[,c(1:3, 370:424, 540)]
#raw intensities
#df_55 <- df[,c(1:3, 255:310, 540)]


colnames(df_55) <- as.character(df_55[3,])
df_55_final <- df_55[-c(1,2,3,4),]

df_55_final <- df_55_final %>%
  rename(Gene.Name = `Gene Name`)



################################################################################
################ cRAP REMOVAL + GENE NAME UPDATE ###############################
################################################################################

## Output directory
#dir.create(here("outputs", "data-processing"), recursive = TRUE, showWarnings = FALSE)

################################################################################
## 1) INPUT DATA
################################################################################

# IMPORTANT:
# df_55_final must already:
# - have correct colnames
# - have "-" replaced by NA
# - have intensity columns numeric
# - contain at least: Accession, Gene Name

proteinGroups <- df_55_final

################################################################################
## 2) cRAP CONTAMINANTS (QC + FILTER)
################################################################################

# Identify cRAP proteins (MaxQuant convention: CON__)
cRAP.proteins <- proteinGroups %>%
  filter(grepl("cRAP", Accession))

# Save cRAP list for QC
write.csv(
  cRAP.proteins,
  here("outputs", "data-processing", "proteinGroups_cRAP.csv"),
  row.names = FALSE
)

write.csv(
  cRAP.proteins,
  here("outputs", "data-processing", "proteinGroups_cRAP_raw.csv"),
  row.names = FALSE
)

# Remove cRAP from main table
d1 <- proteinGroups %>%
  filter(!grepl("cRAP", Accession))

################################################################################
## 3) OPTIONAL: KERATIN REMOVAL (only if Description exists)
################################################################################

if ("Description" %in% colnames(d1)) {
  d1 <- d1 %>%
    filter(!grepl("keratin", Description, ignore.case = TRUE))
}

################################################################################
## 4) MAKE UNIQUE GENE / PROTEIN IDENTIFIERS
################################################################################

# DEP requires unique gene names
d1 <- make_unique(
  d1,
  "Accession",
  "Gene.Name",
  delim = ";"
)

################################################################################
## 5) GENE SYMBOL UPDATE (HGNC)
################################################################################

# Load fixed gene name map (same as paper)
load(here("data", "genenames_update_20220816.RData"))

update_d1 <- checkGeneSymbols(
  d1$Gene.Name,
  species = "human",
  map = genenames_newest,
  unmapped.as.na = FALSE
)

# Manual conflict resolution (paper-consistent)
# predpokladáme, že genenames_newest má stĺpce Symbol a Approved.symbol
gene_map <- setNames(genenames_newest$Approved.Symbol, genenames_newest$Symbol)

# funkcia na výber správneho génu
resolve_gene <- function(x, gene_map) {
  # rozdel podľa oddeľovačov (///, ;, ?)
  parts <- unlist(strsplit(x, "\\s*(///|;|\\?)\\s*"))
  
  # nahraď podľa mapy
  parts_mapped <- sapply(parts, function(g) {
    if(g %in% names(gene_map)) {
      return(gene_map[g])
    } else {
      return(g)
    }
  })
  
  # vyber prvý po mapovaní
  return(parts_mapped[1])
}

# aplikuj na celý stĺpec
update_d1$Suggested.Symbol <- sapply(update_d1$Suggested.Symbol, resolve_gene, gene_map)

update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "C12orf75 /// FSTL1"] <- "C12orf75"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "LNPK /// NUSAP1"] <- "LNPK"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "SARS1 /// SARS2"] <- "SARS1"
update_d1$Suggested.Symbol[update_d1$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"

# Check duplicates (QC)
update_d1$Suggested.Symbol[duplicated(update_d1$Suggested.Symbol)]

# Merge updated symbols back
d1 <- left_join(
  d1,
  update_d1,
  by = c("Gene.Name" = "x")
)

d1$Suggested.Symbol[duplicated(d1$Suggested.Symbol)]

duplicates <- d1 %>% filter(Gene.Name %in% d1$Gene.Name[duplicated(d1$Gene.Name)])
dup_genes <- d1 %>% 
  group_by(Gene.Name) %>% 
  filter(n() > 1) %>% 
  arrange(Gene.Name)

d1 <- d1[!duplicated(d1), ]

colnames(d1)[colnames(d1) == "245_NA"] <- "245_SE"

################################################################################
## 6) FINAL OUTPUT
################################################################################

write.csv(
  d1,
  here("outputs", "02_data-cleaning", "updated_proteinGroups_SEC_only.csv"),
  row.names = FALSE
)

write.csv(
  d1,
  here("outputs", "02_data-cleaning", "updated_proteinGroups_SEC_only_raw.csv"),
  row.names = FALSE
)

################################################################################
## 7) QC SUMMARY
################################################################################

message("Initial proteins: ", nrow(proteinGroups))
message("Removed cRAP proteins: ", nrow(cRAP.proteins))
message("Final proteins after cleaning: ", nrow(d1))

#####################  RNAseq_S2_table ##############################

load(here("data", "RNAseq_S2_table.txt"))
d2 <- RNAseq_S2_table

# Update gene names
update_d2 <- checkGeneSymbols(d2$Cell.type, species = "human", map = genenames_newest, unmapped.as.na = FALSE)

update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "SARS1 /// SARS2"] <- "SARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "GPAT3 /// LPCAT1"] <- "GPAT3"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MTARC1 /// MARCHF1"] <- "MARCHF1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "TAFAZZIN /// WWTR1"] <- "TAFAZZIN"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "CCN3 /// PLXNA1 /// RPL10"] <- "CCN3"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "C11orf98 /// LBHD1"] <- "C11orf98"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "B3GNT2 /// B4GAT1"] <- "B4GAT1"  
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "BHLHE40 /// CENPX"] <- "BHLHE40"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "KAT14 /// PET117"] <- "KAT14"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "IRF4 /// PWWP3A"] <- "PWWP3A"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MTARC2 /// MARCHF2"] <- "MARCHF2"

# Duplicates solving
update_d2$Suggested.Symbol[duplicated(update_d2$Suggested.Symbol)] # check for the duplicates

d2 <- left_join(d2, update_d2, by = c("Cell.type" = "x"))

d2$Suggested.Symbol[duplicated(d2$Suggested.Symbol)] # check for the duplicates
# For these cases delete the older version of the protein
# SRSF10 is there 4x, so deleted
d2 <- d2[ !(d2$Cell.type %in% c("HNRNPU-AS1",  "CASC7",  "LSMD1", "STRA13")), ]
d2 <- d2[ !(d2$Cell.type %in% c("SRSF10")), ]

write.csv(d2, here('outputs', '02_data-cleaning', 'updated_RNAseq_S2_table.csv'))

#####################  Human Cell Map ##############################

d3 <- `preys-latest`

# Update gene names
update_d3 <- checkGeneSymbols(d3$symbol, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d3[update_d3$Approved == FALSE, ])

update_d3$Suggested.Symbol[update_d3$Suggested.Symbol == "SEPTIN2 /// SEPTIN6"] <- "SEPTIN2"
update_d3$Suggested.Symbol[update_d3$Suggested.Symbol == "EPRS1 /// QARS1"] <- "QARS1"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MPHOSPH6 /// PALS2"] <- "PALS2"
update_d2$Suggested.Symbol[update_d2$Suggested.Symbol == "MT-CO2 /// PTGS2"] <- "MT-CO2"

d3 <- left_join(d3, update_d3, by = c("symbol" = "x"))

d3$Suggested.Symbol[duplicated(d3$Suggested.Symbol)] # duplicates

d3 <- d3[ !(d3$symbol %in% c("C10orf12", "KIAA0754")), ]

write.csv(d3, here('outputs', '02_data-cleaning', 'updated_preys-latest.csv'))


#################  MISEV protein categories #########################

d4 <- MISEV2018_protein_categories

# Transform to long-format table
colnames(d4) <- c("cat1", "cat2", "cat3", "cat4", "cat5")
d4 <- pivot_longer(d4, cat1:cat5, names_to = "category", values_to = "marker")
d4 <- d4 %>%
  arrange(category)
d4 <- d4[!d4$marker == "",]

# Update gene names
update_d4 <- checkGeneSymbols(d4$marker, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d4[update_d4$Approved == FALSE, ])
#update_d4$Suggested.Symbol <- sapply(strsplit(update_d4$Suggested.Symbol," /// "), `[`, 1)

d4 <- left_join(d4, update_d4, by = c("marker" = "x"))
d4$Suggested.Symbol[duplicated(d4$Suggested.Symbol)]

d4 <- d4[ !(d4$marker %in% c("IL27", "ITGA2B")), ] # ITGA2B is there 4x, so deleted

write.csv(d4, here('outputs', '02_data-cleaning', 'updated_MISEV_protein_categories.csv'))

#################  Izar (RNAseq) markers #########################

d5 <- Izar_cell_markers_updated

# Update gene names
update_d5 <- checkGeneSymbols(d5$Gene.Name, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d5[update_d5$Approved == FALSE, ])

d5 <- left_join(d5, update_d5, by = c("Gene.Name" = "x"))
d5$Suggested.Symbol[duplicated(d5$Suggested.Symbol)]

write.csv(d5, here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated.csv'))




d4 <- read.delim(
  here("data", "MISEV2018_protein_categories.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Transform to long-format table
colnames(d4) <- c("cat1", "cat2", "cat3", "cat4", "cat5")
d4 <- pivot_longer(d4, cat1:cat5, names_to = "category", values_to = "marker")
d4 <- d4 %>%
  arrange(category)
d4 <- d4[!d4$marker == "",]

# Update gene names
update_d4 <- checkGeneSymbols(d4$marker, species = "human", map = genenames_newest, unmapped.as.na = FALSE)
nrow(update_d4[update_d4$Approved == FALSE, ])
#update_d4$Suggested.Symbol <- sapply(strsplit(update_d4$Suggested.Symbol," /// "), `[`, 1)

d4 <- left_join(d4, update_d4, by = c("marker" = "x"))
d4$Suggested.Symbol[duplicated(d4$Suggested.Symbol)]

d4 <- d4[ !(d4$marker %in% c("IL27", "ITGA2B")), ] # ITGA2B is there 4x, so deleted

write.csv(d4, here('outputs', 'data-processing', 'updated_MISEV_protein_categories.csv'))





################################################################################
######################## MISEV MARKERS MAPPING – SEC-only ######################
################################################################################

library(here)
library(dplyr)
library(readr)
library(stringr)
library(ComplexHeatmap)
library(grid)

# Output dir
dir.create(here("outputs", "MISEV-markers-mapping"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1) INPUTS
################################################################################

# Cleaned proteomics table (output from your cleaning script)
prot_path <- here("outputs", "data-processing", "updated_proteinGroups_SEC_only.csv")

# MISEV categories table (use the one you have; pick ONE of these and delete the other)
# Option A: if you already created it in outputs (from the paper cleaning step)
misev_path <- here("outputs", "data-processing", "updated_MISEV_protein_categories.csv")
# Option B: if you keep it in data/
#misev_path <- here("data", "MISEV2018_protein_categories.txt")

d1 <- read_csv(prot_path, show_col_types = FALSE)
data.MISEV <- read_csv(misev_path, show_col_types = FALSE)

################################################################################
# 2) CREATE A STABLE "UPDATED GENE SYMBOL" COLUMN IN YOUR PROTEOMICS TABLE
################################################################################

# After your HGNC join, d1 should have:
# - Gene.Name (original)
# - Suggested.Symbol (HGNC suggested; can be NA)
# We build one "final" gene symbol used everywhere downstream.

d1 <- d1 %>%
  mutate(Gene.updated = if_else(is.na(Suggested.Symbol) | Suggested.Symbol == "",
                                Gene.Name,
                                Suggested.Symbol))

################################################################################
# 3) IDENTIFY SAMPLE / INTENSITY COLUMNS (SEC-only samples)
################################################################################

# Your sample columns are typically like "041_RP", "056_SE", ...
sample_cols <- names(d1)[grepl("_(SE|RP|RT|NA)$", names(d1))]

# If this returns 0 (naming differs), fallback: keep all numeric-like columns except metadata
if (length(sample_cols) == 0) {
  meta_cols <- c("Accession", "Gene.Name", "Suggested.Symbol", "Gene.updated",
                 "Approved", "ID", "name")
  sample_cols <- setdiff(names(d1), meta_cols)
}

################################################################################
# 4) MAP MISEV MARKERS ONTO YOUR DATA
################################################################################

# data.MISEV is expected to contain:
# - category
# - marker
# - Suggested.Symbol  (updated marker symbol)
# If your MISEV table has different column names, change here accordingly.

data.merged <- data.MISEV %>%
  left_join(
    d1 %>% select(Accession, Gene.Name, Gene.updated, all_of(sample_cols)),
    by = c("marker" = "Gene.updated")
  )

# Export the full mapping table (including markers not found in your data)
write_csv(
  data.merged,
  here("outputs", "MISEV-markers-mapping", "04_MISEV-markers-mapped_LFQlog2_all-markers.csv")
)

################################################################################
# 5) FILTER: keep only markers that were found in your dataset
################################################################################

# Keep only rows where at least one sample has a non-missing value
data.present <- data.merged %>%
  filter(if_any(all_of(sample_cols), ~ !is.na(.x))) %>%
  # coerce sample cols to numeric (important if some columns came as character)
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x))))

# Remove rows that are all zeros or all NA across samples
data.present <- data.present %>%
  rowwise() %>%
  mutate(sum_nonzero = sum(c_across(all_of(sample_cols)) > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(sum_nonzero > 0) %>%
  select(-sum_nonzero)

################################################################################
# 6) BUILD BINARY MATRIX (present/absent)
################################################################################

bin_df <- data.present %>%
  mutate(across(all_of(sample_cols), ~ if_else(!is.na(.x) & .x > 0, 1, 0)))

# Heatmap matrix: rows = markers, cols = patients
mat <- as.matrix(bin_df %>% select(all_of(sample_cols)))
rownames(mat) <- bin_df$marker  # marker symbol from MISEV table

################################################################################
# 7) ANNOTATIONS (SEC-only)
################################################################################

# Row annotation: MISEV category
category <- factor(bin_df$category)
row_ha <- rowAnnotation(category = category)

category <- factor(bin_df$category)
levels(category) <- c("1 - Transmembrane in PM", "2 - Cytosolic",
                      "3 - Non-EV co-isolated", "4 - Transmembrane other",
                      "5 - Secreted")
category_colors <- c(
  "1 - Transmembrane in PM" = "blue",  # tmavozelená
  "2 - Cytosolic" = "#55A868",  # oranžová
  "3 - Non-EV co-isolated" = "yellow",  # fialová
  "4 - Transmembrane other" = "orange",  # ružová
  "5 - Secreted" = "red"   # olivová
)
row_ha <- rowAnnotation(
  category = category,
  col = list(category = category_colors),
  annotation_legend_param = list(
    title = "MISEV protein category"
  )
  )

# Column annotation: phenotype from column name suffix (SE/RP/RT)
#phenotype <- str_match(colnames(mat), "_([A-Z]+)$")[, 4]
#phenotype <- factor(phenotype, levels = c("SE", "RP", "RT", "NA"))
phenotype_colors <- c(
  "SE" = "turquoise3",  # zelená
  "RP" = "mediumvioletred",  # oranžová
  "RT" = "pink",  # 
  "NA" = "grey70"
)

response <- sub(".*_", "", colnames(mat))
response <- factor(response, levels = c("SE", "RP", "RT", "NA"))

ord <- order(response, colnames(mat))
mat_ord <- mat[, ord, drop = FALSE]
phenotype <- response[ord]  # ak používaš premennú phenotype, nech sedí s novým poradím

col_ha <- HeatmapAnnotation(
  phenotype = phenotype,
  col = list(phenotype = phenotype_colors),
  annotation_name_side = "right",
  annotation_legend_param = list(
    title = "Treatment response"
  )
  )

################################################################################
# 8) PLOT HEATMAP
################################################################################

pdf(here("outputs", "MISEV-markers-mapping", "04_MISEV-markers-heatmap_SEC-only.pdf"),
    width = 12, height = 12)

Heatmap(
  mat,
  col = c("0" = "white", "1" = "black"),
  name = "Protein present",
  right_annotation = row_ha,
  top_annotation = col_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_split = category,
  show_row_dend = FALSE,
  row_title = NULL,          # ⬅️ zruší "title" riadkových splitov
  row_names_gp = grid::gpar(fontsize = 0),
  row_title_gp = grid::gpar(fontsize = 0)
)

dev.off()

################################################################################
# 9) EXPORT BINARY MATRIX USED FOR HEATMAP
################################################################################

write.csv(
  mat,
  here("outputs", "MISEV-markers-mapping", "04_MISEV-markers-mapped_binary-matrix_SEC-only.csv"),
  row.names = TRUE
)


library(UpSetR)
library(here)

# mat01: riadky = gény, stĺpce = pacientky, hodnoty 0/1
# napr. mat01 <- data.heatmap2 (ako si mala v predchádzajúcom skripte)

# istota: všetko na numeric 0/1
#mat[] <- lapply(mat, function(x) as.integer(as.character(x)))
#mat[is.na(mat)] <- 0

# UpSetR chce data.frame, kde stĺpce sú množiny a hodnoty 0/1
upset_df <- as.data.frame(mat[, sample_cols, drop = FALSE])
upset_df[, sample_cols][is.na(as.matrix(upset_df[, sample_cols]))] <- 0

# kreslenie – pri 55 stĺpcoch odporúčam obmedziť nsets (napr. top 20)
pdf(here("outputs", "patient-intersection", "upset_SEC_patients_top20.pdf"),
    width = 24, height = 20)
upset(upset_df,
      sets = sample_cols,
      nsets = min(20, length(sample_cols)),
      nintersects = 40,
      order.by = "freq")
dev.off()


gene_presence_n <- rowSums(upset_df[, sample_cols, drop = FALSE])
core75_idx <- gene_presence_n >= 42
upset_df_core75 <- upset_df[core75_idx, sample_cols, drop = FALSE]

dir.create(here("outputs", "patient-intersection"), recursive = TRUE, showWarnings = FALSE)
pdf(here("outputs", "patient-intersection",
         "upset_SEC_core75_patients.pdf"),
    width = 12, height = 10)

upset(
  upset_df_core75,
  nsets = min(20, ncol(upset_df_core75)),
  nintersects = 40,
  order.by = "freq",
  main.bar.color="olivedrab3",
  point.size = 2,
  set_size.show = FALSE,
  text.scale = c(
    2,  # ⬅️ čísla nad hlavnými barmi (intersection size)
    2,  # čísla pri set size baroch
    2,  # názvy setov
    2,  # intersection labels
    1,# x-axis ticks
    2 # y-axis ticks
  )
  )

dev.off()


