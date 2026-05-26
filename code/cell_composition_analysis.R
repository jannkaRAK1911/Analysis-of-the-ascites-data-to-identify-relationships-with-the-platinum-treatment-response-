################################################################################
######################## Cell composition analysis #############################
################################################################################

################################################################################
# LIBRARIES
################################################################################

# Load libraries for data wrangling, visualization, compositional data analysis (CoDA), 
# repeated measures correlation, and multivariate statistics.
library(here)
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(rmcorr)
library(ggprism)
library(compositions)
library(vegan)
library(ggpubr)
library(stringr)
library(tibble)
library(readxl)

# Create output directory for RNA-seq and composition comparisons
dir.create(here('outputs', '5_RNAseq-comparison'), showWarnings = FALSE)

# Load patient ID mapping file
id_map <- read_excel(here("data", "id_map.xlsx"))

# Load cell marker gene expression data
RNAseq_markers <- read.csv(here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated.csv'))

# Filter the dataset to include only specific cell-type marker genes defined by Izar et al.
RNAseq_markers <- RNAseq_markers[RNAseq_markers$Gene.Name %in% c("CLDN3", "FOLR1", "ELF3", "CLDN4", "EPCAM", "TACSTD2", "KRTCAP3", "MMP7", "CLDN7", "SOX17",
                                                                 "PRSS22", "C1orf186", "TSPAN1", "CKB", "S100A14", "LCN2", "SMIM22", "KLF5", "TMEM139",
                                                                 "FXYD3", "PRSS8", "IGF2", "LYNX1", "WNT7A", "SPINT2", "KCNK15", "LYPD1", "SCNN1A", "MEDAG",
                                                                 "FGF7", "COL1A1", "DCN", "SERPINE1", "COL1A2", "VCAM1", "EMILIN1", "NID2", "BDKRB1", "COL3A1",
                                                                 "RGS4", "CDH11", "POSTN", "GPC3", "COL5A1", "CCDC80", "COL5A2", "LOX", "PROCR", "CALB2", "SERPINB2",
                                                                 "TDO2", "FILIP1L", "COL8A1", "PTGER3", "DPP4", "PLS3", "AIF1", "MS4A4A", "MS4A6A", "LY86", "LRRC25",
                                                                 "PILRA", "STAB1", "VSIG4", "CD14", "CD68", "MS4A7", "SLCO2B1", "LILRB2", "CYBB", "ASGR1", "C5AR1",
                                                                 "CD84", "CSF1R", "FCGR1A", "CD163", "ADAP2", "C1QC", "IGSF6", "SLC11A1", "SCIMP", "CD33",
                                                                 "RP11-290F20.3", "FPR1"), ]

# Clean up column names and unify cell type nomenclature (rename CAFs to fibroblast)
RNAseq_markers$X <- NULL
RNAseq_markers$cell.type[RNAseq_markers$cell.type == "CAFs"] <- "fibroblast"


######################################################################
# DATA PREPARATION FOR EV PROTEIN INTENSITIES
######################################################################

# Identify columns containing raw protein intensities ending with specific sample tags
intensity_cols <- names(d1)[grepl("_(SE|RP|RT)$", names(d1))]

# Merge main protein dataset with filtered cell type markers
SEC_filtered <- d1 %>%
  inner_join(RNAseq_markers, by = "Suggested.Symbol")

# Remove missing values and map cell types to cleaner factor levels
SEC_filtered <- na.omit(SEC_filtered) %>% 
  mutate(cell.type = factor(cell.type, levels = c("epithelial_cancer_cells", "fibroblast", "macrophages")))
levels(SEC_filtered$cell.type) <- list("malignant" = "epithelial_cancer_cells", "fibroblast" = "fibroblast", "macrophage" = "macrophages")

# Subset data to retain only ID, cell type, and intensity metrics
SEC_filtered_intensity <- SEC_filtered %>%
  select(ID, cell.type, all_of(intensity_cols))

# Set gene IDs as row names
rownames(SEC_filtered_intensity) <- SEC_filtered_intensity$ID
SEC_filtered_intensity$ID <- NULL

# Reshape data into long format for ggplot consumption
SEC_long_intensity <- SEC_filtered_intensity %>%
  rownames_to_column(var = "gene") %>%  
  pivot_longer(
    cols = all_of(intensity_cols),
    names_to = "sample",
    values_to = "intensity"
  )

# Filter out missing records represented as dashes and cast intensity to numeric values
SEC_long_intensity <- SEC_long_intensity %>%
  filter(intensity != "-") %>% 
  mutate(intensity = as.numeric(intensity))

# Extract original patient numeric identifiers from sample string and join map to get updated IDs
SEC_long_intensity <- SEC_long_intensity %>%
  mutate(
    old_id = str_extract(sample, "^[^_]+"), 
    old_id = as.numeric(old_id)             
  ) %>%
  left_join(id_map, by = c("old_id" = "old_id"))


######################################################################
# VISUALIZATION: EV-BASED COMPOSITION
######################################################################

# Plot 1: Total relative cell type composition across all individual patients based on EV data
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_rawIntensities.pdf'), width = 38, height = 15)
SEC_long_intensity %>%
  group_by(new_id, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup() %>%
  mutate(sample_factor = factor(new_id, levels = unique(new_id))) %>%
  ggplot(aes(x = sample_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Relative cell type composition across patients",
    subtitle = "Based on EV proteins",
    x = "Sample",
    y = "Relative Intensity",
    fill = "Cell Type"
  ) +
  scale_fill_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  )) +
  theme_minimal(base_size = 18) +
  theme( 
    plot.title = element_text(size = 70, face = "bold"),
    plot.subtitle = element_text(size=45, face = "bold"),
    axis.title = element_text(size = 40), 
    axis.text = element_text(size = 38), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.title = element_text(size = 40), 
    legend.text = element_text(size = 38),
    legend.position = "right" 
  )
dev.off()

# Plot 2: Relative cell type composition broken down by clinical treatment response groups (RP, RT vs SE)
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_relativeIntensities_response.pdf'), width = 30, height = 15)

# Structure and calculate relative fractions by response categories
df_SEC <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  mutate(patient_id = as.numeric(patient_id)) %>%
  left_join(id_map, by = c("patient_id" = "old_id")) %>%
  filter(!is.na(intensity)) %>%
  group_by(new_id.x, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(new_id.x, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup() %>%
  mutate(patient_factor = factor(new_id.x, levels = unique(new_id.x)))

# Define shared plot styling parameters
common_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold", hjust = 0),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 38),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 38),
    strip.text = element_text(size = 40, face = "bold")
  )

# Upper panels: Stacked barcharts for partial response (RP) and stable treatment (RT) cohorts
p_top <- df_SEC %>%
  filter(response %in% c("RP", "RT")) %>%
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~response, scales = "free_x") +
  labs(x = "Patient", y = "Relative intensity", fill = "Cell type") +
  scale_fill_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  )) +
  common_theme

# Lower panel: Stacked barchart spanning the full layout for side-effects/other group (SE)
p_bottom <- df_SEC %>%
  filter(response == "SE") %>%
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  labs(x = "Patient", y = "Relative intensity", fill = "Cell type") +
  ggtitle("SE") +
  common_theme +
  scale_fill_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  )) +
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))

# Assemble top and bottom tracks using patchwork architecture
final_plot <- p_top / p_bottom +
  plot_layout(guides = "collect", heights = c(1, 1)) & 
  theme(legend.position = "right")

# Output plot with global master annotations
final_plot +
  plot_annotation(
    title = "Relative cell type composition across patients by treatment response",
    subtitle = "Based on EV proteins",
    theme = theme(
      plot.title = element_text(size = 60, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 40, face = "bold")
    )
  )
dev.off()


###################################################################
# STATISTICAL ANALYSIS & CODA TRANSFORMATIONS (EV)
###################################################################

# Consolidate dynamic proportions for statistical testing
df_test <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  group_by(patient_id, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(patient_id, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup() %>%
  mutate(patient_id = as.numeric(patient_id))

# Define tiny imputation constant to resolve zero-values in compositional signatures
pseudocount <- 1e-6

# Apply Centered Log-Ratio (CLR) transformation to handle compositional constraints
df_test_clr <- df_test %>%
  left_join(
    id_map %>% select(old_id, new_id),
    by = c("patient_id" = "old_id")
  ) %>%
  group_by(new_id) %>%
  mutate(intensity.relative = intensity.relative + pseudocount) %>%
  mutate(clr_value = log(intensity.relative / exp(mean(log(intensity.relative))))) %>%
  ungroup()

# Plot boxplots comparing CLR-transformed values with pairwise Wilcoxon tests
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_boxplot-SEC_relativeIntensities_CLR.pdf'), width = 10, height = 6)

df_test_clr %>%
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  ggplot(aes(x = response, y = clr_value, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~cell.type, scales = "free_y") +
  stat_compare_means(
    comparisons = list(c("RP", "RT"), c("RP", "SE"), c("RT", "SE")), # added dynamic list for consistency
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE
  ) +
  scale_fill_manual(values = c("RT" = "#D55E00",
                               "RP" = "#CC79A7",
                               "SE" = "#0072B2")) +
  scale_color_manual(values = c("RT" = "#D55E00",
                                "RP" = "#CC79A7",
                                "SE" = "#0072B2")) +
  labs(
    title = "Comparison of CLR-transformed cell type composition",
    subtitle = "Estimated from EV protein intensities",
    x = "Response",
    y = "CLR transformed relative proportion",
    fill = "Response",
    color = "Response"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 16, hjust = 0, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
dev.off()


###########################################################################
# EXTENDED CELL TYPES ANALYSIS (FLOW CYTOMETRY DATA) 
###########################################################################

# Load gold-standard profiling information acquired via Flow Cytometry (FC)
df_all_cells <- read_excel(here("data", "flowcytometry_data.xlsx"))
df_all_cells <- df_all_cells %>%
  rename(patient_id = `Patient ID`) %>%
  left_join(id_map, by = c("patient_id" = "old_id"))

# Plot 3: Expanded cell populations overview across individual patient landscapes
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_rawIntensities_all_cells.pdf'), width = 38, height = 15)

df_all_cells %>%
  pivot_longer(
    cols = c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
             "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells"),
    names_to = "cell.type",
    values_to = "proportion"
  ) %>%
  mutate(proportion = proportion / 100) %>% # Normalise scale from percentages to fractions
  mutate(patient_factor = factor(new_id, levels = sort(unique(new_id)))) %>%
  ggplot(aes(x = patient_factor, y = proportion, fill = cell.type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Malignant cells" = "#F8766D",
    "Fibroblasts" = "#00BA38",
    "Macrophages" = "#619CFF",
    "Monocytes" = "#3F4384",
    "B cells" = "#C79CC5",
    "CD4 T cells" = "#A45DA8",
    "CD8 T cells" = "#823F9B",
    "Neutrophils" = "#FAB038",
    "NK cells" = "#999999",
    "Other cells" = "black"
  )) +
  labs(
    title = "Relative cell type composition across patients",
    subtitle = "Based on FC data",
    x = "Patient",
    y = "Relative proportion",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 18) +
  theme( 
    plot.title = element_text(size = 70, face = "bold"),
    plot.subtitle = element_text(size=45, face = "bold"),
    axis.title = element_text(size = 40), 
    axis.text = element_text(size = 38), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.title = element_text(size = 40), 
    legend.text = element_text(size = 38),
    legend.position = "right"
  )
dev.off()

# Plot 4: Flow Cytometry comprehensive populations stratified by therapeutic group split
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_relativeIntensities_response_all_cells.pdf'), width = 30, height = 15)

df_long <- df_all_cells %>%
  pivot_longer(
    cols = c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
             "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells"),
    names_to = "cell.type",
    values_to = "proportion"
  ) %>%
  mutate(proportion = proportion / 100) %>%
  mutate(patient_factor = factor(new_id, levels = sort(unique(new_id))))

# Construct composite graphs mapping responses for flow cytometry metrics
p_top <- df_long %>%
  filter(response %in% c("RP", "RT")) %>%
  ggplot(aes(x = patient_factor, y = proportion, fill = cell.type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~response, scales = "free_x") +
  scale_fill_manual(values = c(
    "Malignant cells" = "#F8766D",
    "Fibroblasts" = "#00BA38",
    "Macrophages" = "#619CFF",
    "Monocytes" = "#3F4384",
    "B cells" = "#C79CC5",
    "CD4 T cells" = "#A45DA8",
    "CD8 T cells" = "#823F9B",
    "Neutrophils" = "#FAB038",
    "NK cells" = "#999999",
    "Other cells" = "black"
  )) +
  labs(x = "Patient", y = "Relative proportion", fill = "Cell Type") +
  common_theme

p_bottom <- df_long %>%
  filter(response == "SE") %>%
  ggplot(aes(x = patient_factor, y = proportion, fill = cell.type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Malignant cells" = "#F8766D",
    "Fibroblasts" = "#00BA38",
    "Macrophages" = "#619CFF",
    "Monocytes" = "#3F4384",
    "B cells" = "#C79CC5",
    "CD4 T cells" = "#A45DA8",
    "CD8 T cells" = "#823F9B",
    "Neutrophils" = "#FAB038",
    "NK cells" = "#999999",
    "Other cells" = "black"
  )) +
  labs(x = "Patient", y = "Relative proportion", fill = "Cell Type") +
  ggtitle("SE") +
  common_theme +
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5))

final_plot_fc <- p_top / p_bottom +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")

final_plot_fc +
  plot_annotation(
    title = "Relative cell type composition across patients by response",
    subtitle = "Based on FC data",
    theme = theme(
      plot.title = element_text(size = 55, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 40, face = "bold")
    )
  )
dev.off()


###################################################################
# STATISTICAL ANALYSIS & TRANSFORMATIONS (FC)
###################################################################

# Apply compositional centered log-ratio (CLR) conversion to full flow cytometry arrays
df_clr <- df_long %>%
  group_by(new_id) %>%   
  mutate(proportion = proportion + pseudocount) %>%
  mutate(clr_value = log(proportion / exp(mean(log(proportion))))) %>%
  ungroup()

# Set testing comparison matrices
comparisons <- list(c("RP", "RT"), c("RP", "SE"), c("RT", "SE"))

# Plot 5: Matrix boxplot grids of all 10 cell subsets using CLR signals
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_boxplot-SEC_CLR_all_cells.pdf'), width = 15, height = 12)

df_clr %>%
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  ggplot(aes(x = response, y = clr_value, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~cell.type, scales = "free_y", nrow = 2, ncol = 5) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE,
    size = 5
  ) +
  scale_fill_manual(values = c("RT" = "#D55E00",
                               "RP" = "#CC79A7",
                               "SE" = "#0072B2")) +
  
  scale_color_manual(values = c("RT" = "#D55E00",
                                "RP" = "#CC79A7",
                                "SE" = "#0072B2")) +
  labs(
    title = "Comparison of CLR-transformed cell type composition",
    subtitle = "Based on the FC data",
    x = "Response",
    y = "CLR transformed relative proportion",
    fill = "Response",
    color = "Response"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10),
    panel.spacing = unit(2, "lines"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 23, hjust = 0, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
dev.off()


###################################################################
# REPEATED MEASURES CORRELATION: FLOW CYTOMETRY VS EVS 
###################################################################

# Isolate matching major categories (Malignant, Fibroblast, Macrophage) from flow dataset
df_cells_3 <- df_all_cells %>%
  select(new_id, response, "Malignant cells", Fibroblasts, Macrophages) %>%
  rename(
    malignant = "Malignant cells",
    fibroblast = Fibroblasts,
    macrophage = Macrophages,
    patient_id = new_id
  )

# Normalize selected sub-composition to re-establish sum-to-100% constant closure
df_cells_3_norm <- df_cells_3 %>%
  mutate(
    total = malignant + fibroblast + macrophage,
    malignant = malignant / total * 100,
    fibroblast = fibroblast / total * 100,
    macrophage = macrophage / total * 100
  ) %>%
  select(-total)

# Melt scaled flow datasets to long layout variations
df_cells_3_long <- df_cells_3_norm %>%
  pivot_longer(
    cols = c(macrophage, fibroblast, malignant),
    names_to = "cell.type",
    values_to = "FC"
  ) %>%
  arrange(patient_id, response, cell.type)

# Format Extracellular Vesicle (EV) protein relative subsets for direct pair pairing
df_test_percentages <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  group_by(new_id, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(new_id, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup()

df_test_percentages_long <- df_test_percentages %>%
  select(new_id, response, cell.type, intensity.relative) %>%
  rename(patient_id = new_id, EV = intensity.relative)

# Filter cohort identifiers present in both diagnostic tracks
common_ids <- intersect(df_cells_3_long$patient_id, df_test_percentages_long$patient_id)

df1 <- df_cells_3_long %>% filter(patient_id %in% common_ids)
df2 <- df_test_percentages_long %>% filter(patient_id %in% common_ids)

# Combine datasets by patient ID and cell category to construct rmcorr evaluation structures
df_rmcorr <- left_join(df2, df1, by = c("patient_id", "cell.type"))

# Compute global Repeated Measures Correlation coefficient evaluating tracking alignment 
my.rmc.cells.EV <- rmcorr(
  participant = cell.type,
  measure1 = FC,
  measure2 = EV,
  dataset = df_rmcorr
)

# Extract baseline Spearman rank variables to construct descriptive visualization labels
corr_df <- df_rmcorr %>%
  group_by(cell.type) %>%
  summarise(
    r = cor(FC, EV * 100, method = "spearman", use = "complete.obs"),
    p = cor.test(FC, EV * 100)$p.value
  ) %>%
  mutate(
    r_label = case_when(
      cell.type == "fibroblast" ~ paste0("italic(r)[f] == ", round(r, 2)),
      cell.type == "malignant"  ~ paste0("italic(r)[m] == ", round(r, 2)),
      cell.type == "macrophage" ~ paste0("italic(r)[mac] == ", round(r, 2))
    )
  )

# Plot 6: Combined regression analysis tracking alignment between physical cell metrics and EV profiles
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_rmcorr_FCvsEV_new.pdf'), width = 10, height = 8)
ggplot(df_rmcorr, aes(x = FC, y = EV*100, color = cell.type)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1) +
  geom_text(
    data = corr_df,
    aes(x = 80, y = c(95, 88, 81), label = r_label, color = cell.type),
    parse = TRUE, size = 5, show.legend = FALSE
  ) +
  theme_prism() +
  xlim(0, 100) + ylim(0, 100) +
  labs(x = "FC (%)", y = "EV (%)", title = "Repeated measure correlation") +
  scale_color_manual(values = c("malignant" = "#F8766D",
                                "fibroblast" = "#00BA38",
                                "macrophage" = "#619CFF")) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  ) +
  annotate("text", x = 60, y = 14, label = paste("r =", round(my.rmc.cells.EV$r, digits = 3)), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells.EV$CI[1], digits = 3), "-", round(my.rmc.cells.EV$CI[2], digits = 3)), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 2, label = paste("p =", format(my.rmc.cells.EV$p, scientific = TRUE, digits = 3)), hjust = 0, size = 5)
dev.off()


###################################################################
# MULTIVARIATE STATISTICS: PERMANOVA & PCA BLOCKS 
###################################################################

### 1. Analysis using All Flow Cytometry Cell Types (10 dimensions)
cell_data <- df_all_cells[, c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
                              "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells")]

# Impute zeros, transform via CLR, and construct Aitchison distance metrics
cell_data[cell_data == 0] <- pseudocount
clr_data <- clr(cell_data)
dist_matrix <- dist(clr_data)

# Run PERMANOVA to test structural differences driven by clinical response variations
adonis2(dist_matrix ~ response, data = df_all_cells)

# Compute Principal Component Analysis (PCA) space coordinates
pca <- prcomp(clr_data)
pca_df <- data.frame(pca$x, response = df_all_cells$response)

# Save visual configuration of full FC profiles spatial distribution
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_all_cells_raw.pdf'), width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = response)) +
  scale_color_manual(values = c("RT" = "#D55E00",
                                "RP" = "#CC79A7",
                                "SE" = "#0072B2")) +
  geom_point(size = 6) +
  theme(
    axis.title = element_text(size = 35),   
    axis.text = element_text(size = 29),    
    legend.title = element_text(size = 31), 
    legend.text = element_text(size = 29)   
  )
dev.off()


### 2. Analysis using Major 3 Flow Cytometry Cell Types Sub-composition
cell_data.3 <- df_cells_3_norm[, c("malignant", "fibroblast", "macrophage")]

# Clean composition constraints, convert to Euclidean targets, and test groupings
cell_data.3[cell_data.3 == 0] <- pseudocount
clr_data.3 <- clr(cell_data.3)
dist_matrix.3 <- dist(clr_data.3)

adonis2(dist_matrix.3 ~ response, data = df_cells_3_norm)

# Map dynamic projections using PCA models
pca.3 <- prcomp(clr_data.3)
pca_df.3 <- data.frame(pca.3$x, response = df_cells_3_norm$response)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_cells_3_raw.pdf'), width = 10, height = 8)
ggplot(pca_df.3, aes(x = PC1, y = PC2, color = response)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("RT" = "#D55E00",
                                "RP" = "#CC79A7",
                                "SE" = "#0072B2")) +
  theme(
    axis.title = element_text(size = 35),   
    axis.text = element_text(size = 29),    
    legend.title = element_text(size = 31), 
    legend.text = element_text(size = 29)   
  )
dev.off()


### 3. Analysis using Extracellular Vesicle (EV) Relative Fractions
EV_data <- df_test_percentages_long %>%
  pivot_wider(names_from = cell.type, values_from = EV)

EV_data_cells <- EV_data[, c("malignant", "fibroblast", "macrophage")]
EV_data_cells[EV_data_cells == 0] <- pseudocount

# Perform CLR conversion and evaluate shifts inside matrix networks
clr_data_EV <- clr(EV_data_cells)
dist_matrix_EV <- dist(clr_data_EV)

adonis2(dist_matrix_EV ~ response, data = EV_data)

# Build PCA layout mirroring mapped variance parameters
pca_EV <- prcomp(clr_data_EV)
pca_df_EV <- data.frame(pca_EV$x, response = EV_data$response)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_EVs_raw.pdf'), width = 10, height = 8)
ggplot(pca_df_EV, aes(x = PC1, y = PC2, color = response)) +
  geom_point(size = 6) +
  scale_color_manual(values = c("RT" = "#D55E00",
                                "RP" = "#CC79A7",
                                "SE" = "#0072B2")) +
  theme(
    axis.title = element_text(size = 35),   
    axis.text = element_text(size = 29),    
    legend.title = element_text(size = 31), 
    legend.text = element_text(size = 29)   
  )

dev.off()
devdev.off()
