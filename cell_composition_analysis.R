install.packages("betareg")
install.packages("emmeans")
install.packages("rmcorr")
install.packages("ggprism")
install.packages("lme4")
install.packages("compositions")
install.packages("tidyverse")
library(ComplexHeatmap)
library(here)
library(dplyr)
library(patchwork)
library(ggplot2)
library(betareg)
library(emmeans)
library(tidyr)
library(rmcorr)
library(ggprism)
library(lme4)
library(compositions)
library(vegan)
library(ggpubr)
library(stringr)
library(tibble)
library(tidyverse)
library(readxl)


dir.create(here('outputs', '5_RNAseq-comparison'))

# Load the data:
RNAseq_markers <- read.csv(here('outputs', '02_data-cleaning', 'updated_Izar_cell_markers_updated.csv'))

# Filter just the subset of proteins defined by Izar:
RNAseq_markers <- RNAseq_markers[RNAseq_markers$Gene.Name %in% c("CLDN3", "FOLR1", "ELF3", "CLDN4", "EPCAM", "TACSTD2", "KRTCAP3", "MMP7", "CLDN7", "SOX17",
                                                                 "PRSS22", "C1orf186", "TSPAN1", "CKB", "S100A14", "LCN2", "SMIM22", "KLF5", "TMEM139",
                                                                 "FXYD3", "PRSS8", "IGF2", "LYNX1", "WNT7A", "SPINT2", "KCNK15", "LYPD1", "SCNN1A", "MEDAG",
                                                                 "FGF7", "COL1A1", "DCN", "SERPINE1", "COL1A2", "VCAM1", "EMILIN1", "NID2", "BDKRB1", "COL3A1",
                                                                 "RGS4", "CDH11", "POSTN", "GPC3", "COL5A1", "CCDC80", "COL5A2", "LOX", "PROCR", "CALB2", "SERPINB2",
                                                                 "TDO2", "FILIP1L", "COL8A1", "PTGER3", "DPP4", "PLS3", "AIF1", "MS4A4A", "MS4A6A", "LY86", "LRRC25",
                                                                 "PILRA", "STAB1", "VSIG4", "CD14", "CD68", "MS4A7", "SLCO2B1", "LILRB2", "CYBB", "ASGR1", "C5AR1",
                                                                 "CD84", "CSF1R", "FCGR1A", "CD163", "ADAP2", "C1QC", "IGSF6", "SLC11A1", "SCIMP", "CD33",
                                                                 "RP11-290F20.3", "FPR1"), ]

RNAseq_markers$X <- NULL
RNAseq_markers$cell.type[RNAseq_markers$cell.type == "CAFs"] <- "fibroblast"


patients <- d1$Suggested.Symbol
patients <- data.frame(value = patients)

data.filtered <- read.csv(here('outputs', '03_data-processed', '03_data-processed.csv'))
data.filtered.binary <- data.filtered
data.filtered.binary[is.na(data.filtered.binary)] <- 0

# Map marker data onto our protein data
data.merged <- left_join(RNAseq_markers, d1, by = c("Suggested.Symbol" = "Suggested.Symbol"))

#intensity_cols <- grep("^X", colnames(data), value = TRUE)


# Keep only rows where at least one sample has a non-missing value
data.present <- data.merged %>%
  filter(if_any(all_of(intensity_cols), ~ !is.na(.x))) %>%
  # coerce sample cols to numeric (important if some columns came as character)
  mutate(across(all_of(intensity_cols), ~ suppressWarnings(as.numeric(.x))))

# Remove rows that are all zeros or all NA across samples
data.present <- data.present %>%
  rowwise() %>%
  mutate(sum_nonzero = sum(c_across(all_of(intensity_cols)) > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(sum_nonzero > 0) %>%
  select(-sum_nonzero)

data.present[is.na(data.present)] <- 0

data.merged.binary <- data.present %>%
  mutate(across(all_of(intensity_cols), ~ if_else(!is.na(.x) & .x > 0, 1, 0)))
# všetko >0 = 1
data.merged.binary[, intensity_cols][data.merged.binary[, intensity_cols] > 0] <- 1

# odstrániť proteíny, kde sú všetky 0
data.merged.binary <- data.merged.binary[rowSums(data.merged.binary[, intensity_cols]) != 0, ]

data.heatmap <- as.matrix(data.merged.binary[, intensity_cols])

rownames(data.heatmap) <- data.merged.binary$Suggested.Symbol

RNAseq_subset_binary <- left_join(patients, RNAseq_markers, by = c("value" = "Suggested.Symbol"))
RNAseq_subset_binary <- na.omit(RNAseq_subset_binary) %>% 
  mutate(cell.type = factor(cell.type, levels = c("epithelial_cancer_cells", "fibroblast", "macrophages"))) 
levels(RNAseq_subset_binary$cell.type) <- list("malignant" = "epithelial_cancer_cells", "fibroblast" = "fibroblast", "macrophage" = "macrophages")

RNAseq_subset_binary <- left_join(RNAseq_subset_binary, data.merged.binary, by = c("value" = "Suggested.Symbol.x"))

RNAseq_subset_binary_SEC_mat <- RNAseq_subset_binary %>%
  select(all_of(intensity_cols)) %>%
  as.matrix()
rownames(RNAseq_subset_binary_SEC_mat) <-RNAseq_subset_binary$value

category <- RNAseq_subset_binary$cell.type.x

row_ha <- rowAnnotation(category = category)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison-heatmap-binary.pdf'))
Heatmap(RNAseq_subset_binary_SEC_mat,  col = c("white", "black"),
        right_annotation = row_ha,
        split = RNAseq_subset_binary$cell.type.x,
        cluster_row_slices = FALSE,
        row_title = NULL,
        row_names_gp = gpar(fontsize = 5))
dev.off()


# SEC barplots (upravené podľa intensity_cols)
SEC1 <- RNAseq_subset_binary_SEC_mat <- RNAseq_subset_binary %>%
  select(all_of(intensity_cols)) %>% 
  as.matrix()
rownames(SEC1) <- RNAseq_subset_binary$value

# Vyber len riadky z pôvodných dát, ktoré sú v heatmap matrix
SEC2 <- data.present[data.present$Suggested.Symbol %in% rownames(SEC1), ]

# Odstrániť nepotrebné stĺpce
SEC2 <- SEC2 %>% select(-all_of(setdiff(names(SEC2), c("Suggested.Symbol", intensity_cols))))
rownames(SEC2) <- SEC2$Suggested.Symbol
SEC2$Suggested.Symbol <- NULL

# Zoradiť podľa názvov riadkov
SEC1 <- SEC1[order(row.names(SEC1)), ]
SEC2 <- SEC2[order(row.names(SEC2)), ]

# Naplniť matice hodnotami intenzít podľa binárnej matice
SEC_final <- SEC1
for (i in 1:nrow(SEC_final)) {
  for (j in 1:ncol(SEC_final)){
    if (SEC_final[i,j] == 1) {
      SEC_final[i,j] <- SEC2[i,j]
    }
  }
}

# Pivot do long formátu pre ggplot
SEC_long <- as.data.frame(SEC_final) %>%
  tibble::rownames_to_column("rowname") %>%
  left_join(RNAseq_subset_binary %>% select(value, cell.type.x), 
            by = c("rowname" = "value")) %>%
  tidyr::pivot_longer(cols = all_of(intensity_cols), names_to = 'patient', values_to = 'intensity')

# Prepis názvov pacientov, ak je potrebné
SEC_long$patient <- gsub(pattern = "S", replacement = "patient_", SEC_long$patient)

# Vykreslenie barplotu
library(ggplot2)
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC.pdf'))
SEC_long %>%
  group_by(patient, cell.type.x) %>%
  mutate(mean.intensity = mean(intensity)) %>%
  ungroup() %>%
  group_by(patient) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  mutate(patients_factor = factor(patient, levels = unique(patient))) %>%
  ggplot(aes(x = patients_factor, y = intensity.relative, fill = cell.type.x)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "RNAseq")
dev.off()


######################################################################
########### FINAL SCRIPT #############################################
######################################################################

# colnames(d1)[colnames(d1) == "Suggested.Symbol.x"] <- "Suggested.Symbol"
intensity_cols <- names(d1)[grepl("_(SE|RP|RT)$", names(d1))]
SEC_filtered <- d1 %>%
  inner_join(RNAseq_markers, by = "Suggested.Symbol")

SEC_filtered <- na.omit(SEC_filtered) %>% 
  mutate(cell.type = factor(cell.type, levels = c("epithelial_cancer_cells", "fibroblast", "macrophages")))
levels(SEC_filtered$cell.type) <- list("malignant" = "epithelial_cancer_cells", "fibroblast" = "fibroblast", "macrophage" = "macrophages")

SEC_filtered_intensity <- SEC_filtered %>%
  select(ID, cell.type, all_of(intensity_cols))

rownames(SEC_filtered_intensity) <- SEC_filtered_intensity$ID
SEC_filtered_intensity$ID <- NULL

id_map <- read_excel(here("data","id_map.xlsx"))

SEC_long_intensity <- SEC_filtered_intensity %>%
  rownames_to_column(var = "gene") %>%   # ← tu si uložíš rownames
  pivot_longer(
    cols = all_of(intensity_cols),
    names_to = "sample",
    values_to = "intensity"
  )
SEC_long_intensity <- SEC_long_intensity %>%
  filter(intensity != "-") %>%   # odstráni tie pomlčky
  mutate(intensity = as.numeric(intensity))

SEC_long_intensity <- SEC_long_intensity %>%
  mutate(
    old_id = str_extract(sample, "^[^_]+"), # vytiahne 041
    old_id = as.numeric(old_id)              # zmení na 41
  ) %>%
  left_join(id_map, by = c("old_id" = "old_id"))

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_rawIntensities.pdf'), width = 38, height = 15)
SEC_long_intensity %>%
  group_by(new_id, cell.type) %>%
  mutate(mean.intensity = mean(intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(new_id) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
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
  ))+
  theme_minimal(base_size = 18) +
  theme( 
    plot.title = element_text(size = 70, face = "bold"),
    plot.subtitle = element_text(size=45, face = "bold"),
    axis.title = element_text(size = 40), 
    axis.text = element_text(size = 38), # čísla + sample_factor 
    axis.text.x = element_text(angle = 45, hjust = 1), # otočenie + veľkosť 
    legend.title = element_text(size = 40), # legenda názov 
    legend.text = element_text(size = 38),
    legend.position = "right" 
  )
dev.off()


pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_normalizedIntensities_response.pdf'), width = 25, height = 15)
SEC_long_intensity %>%
  # rozdelenie sample na patient_id a response
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  
  # odstránenie NA a zabezpečenie numeric hodnoty
  filter(!is.na(intensity)) %>%
  
  # výpočet priemernej intenzity pre cell type v rámci pacienta
  group_by(patient_id, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  
  # prepočet na relatívne zastúpenie
  group_by(patient_id, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  
  # faktor pre správne zoradenie pacientov
  mutate(patient_factor = factor(patient_id, levels = unique(patient_id))) %>%
  
  # plot
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~response, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Proteomics by Treatment Response",
    x = "Patient",
    y = "Relative intensity"
  )

SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  filter(!is.na(intensity)) %>%
  
  group_by(patient_id, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  
  group_by(patient_id, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  
  mutate(patient_factor = factor(patient_id, levels = unique(patient_id))) %>%
  
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~response, scales = "free_x", ncol = 2) +
  
  labs(
    title = "Proteomics by Treatment Response",
    x = "Patient",
    y = "Relative intensity",
    fill = "Cell type"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 60, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 38),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 38),
    strip.text = element_text(size = 40, face = "bold")
  )


pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_relativeIntensities_response.pdf'), width = 30, height = 15)
# 📊 príprava dát
df_SEC <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  mutate(
    patient_id = as.numeric(patient_id)
  ) %>%
  left_join(id_map, by = c("patient_id" = "old_id")) %>%
  filter(!is.na(intensity)) %>%
  group_by(new_id.x, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(new_id.x, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup() %>%
  mutate(patient_factor = factor(new_id.x, levels = unique(new_id.x)))


# 🎨 spoločný theme
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

# 🔝 horné panely (RP, RT)
p_top <- df_SEC %>%
  filter(response %in% c("RP", "RT")) %>%
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~response, scales = "free_x") +
  labs(
    x = "Patient",
    y = "Relative intensity",
    fill = "Cell type"
  ) +
  scale_fill_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  ))+
  common_theme

# 🔻 dolný panel (SE cez celú šírku)
p_bottom <- df_SEC %>%
  filter(response == "SE") %>%
  ggplot(aes(x = patient_factor, y = intensity.relative, fill = cell.type)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Patient",
    y = "Relative intensity",
    fill = "Cell type"
  ) +
  ggtitle("SE") +
  common_theme +
  scale_fill_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  ))+
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
  )

# 🧩 spojenie grafov
final_plot <- p_top / p_bottom +
  plot_layout(guides = "collect", heights = c(1, 1)) &  # jedna legenda
  theme(legend.position = "right")

# ➕ pridanie titulku
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


### Statistical analysis

# priprava dát
df_test <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  group_by(patient_id, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(patient_id, response) %>%
  mutate(intensity.relative = mean.intensity / sum(mean.intensity)) %>%
  ungroup()

# test pre každý cell type
results <- df_test %>%
  group_by(cell.type) %>%
  summarise(
    p_value = kruskal.test(intensity.relative ~ response)$p.value
  )

results


library(vegan)
# wide formát: cell types ako stĺpce
df_wide <- df_test %>%
  select(patient_id, response, cell.type, intensity.relative) %>%
  pivot_wider(names_from = cell.type, values_from = intensity.relative)

# distance matrix
dist_mat <- vegdist(df_wide[,-c(1,2)], method = "euclidean")

# PERMANOVA
adonis2(dist_mat ~ response, data = df_wide)

library(ggplot2)
library(ggpubr)
ggboxplot(df_test, x = "response", y = "cell.type") +
  stat_compare_means(method = "kruskal.test")

ggplot(df_test, aes(x = response, y = cell.type)) +
  geom_boxplot() +
  theme_minimal()

comparisons <- list(
  c("RP", "SE"),
  c("RP", "RT"),
  c("SE", "RT")
)

df_test %>%
  filter(cell.type == "malignant") %>%
  ggplot(aes(x = response, y = intensity.relative)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     p.adjust.method = "BH") +
  theme_minimal()


###################################################################
##################### Wilcoxon + boxploty #########################
############################# CLR #################################
###################################################################

df_test <- df_test %>%
  mutate(patient_id = as.numeric(patient_id))

df_test_clr <- df_test %>%
  
  left_join(
    id_map %>% select(old_id, new_id),
    by = c("patient_id" = "old_id")
  ) %>%
  
  # 🔄 CLR transformácia v rámci každej vzorky
  group_by(new_id) %>%
  
  mutate(
    intensity.relative = intensity.relative + pseudocount
  ) %>%
  
  mutate(
    clr_value = log(
      intensity.relative /
        exp(mean(log(intensity.relative)))
    )
  ) %>%
  
  ungroup()

# 📄 PDF
pdf(
  here(
    'outputs',
    '5_RNAseq-comparison',
    '5_RNAseq-comparison_boxplot-SEC_relativeIntensities_CLR.pdf'
  ),
  width = 10,
  height = 6
)

# 📊 plot
df_test_clr %>%
  
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  
  ggplot(aes(x = response, y = clr_value, fill = response)) +
  
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  
  geom_jitter(
    aes(color = response),
    width = 0.1,
    alpha = 0.5,
    size = 1.5
  ) +
  
  facet_wrap(~cell.type, scales = "free_y") +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE
  ) +
  
  scale_fill_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  )) +
  
  scale_color_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  )) +
  
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

df_test_clr %>%
  group_by(cell.type) %>%
  summarise(
    p_value = shapiro.test(clr_value)$p.value
  )
ggplot(df_test_clr, aes(sample = clr_value)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~cell.type, scales = "free")

ggplot(df_test_clr, aes(x = clr_value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~cell.type, scales = "free")

###################################################################
##################### Wilcoxon + boxploty #########################
#################### Non transformed data #########################
###################################################################

df_test %>%
  filter(cell.type == "malignant") %>%
  ggplot(aes(x = response, y = intensity.relative, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Malignant cell proportion across response groups",
    x = "Treatment response",
    y = "Relative intensity",
    fill = "Response",
    color = "Response"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

df_test %>%
  filter(cell.type == "malignant") %>%
  ggplot(aes(x = response, y = intensity.relative, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Malignant cell proportion across response groups",
    x = "Treatment response",
    y = "Relative intensity",
    fill = "Response",
    color = "Response"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

df_test %>%
  filter(cell.type == "malignant") %>%
  ggplot(aes(x = response, y = intensity.relative, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.6, size = 2) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     p.adjust.method = "BH",
                     hide.ns = TRUE) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Malignant cell proportion across response groups",
    x = "Treatment response",
    y = "Relative intensity",
    fill = "Response",
    color = "Response"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )


pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_boxplot-SEC_meanIntensities.pdf'), width = 10, height = 6)
df_test %>%
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  ggplot(aes(x = response, y = mean.intensity, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~cell.type, scales = "free_y") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Comparison of cell type composition across response groups",
    subtitle = "Mean intensity of EV proteins",
    x = "Response",
    y = "Mean Intensity",
    fill = "Response",
    color = "Response"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = -0.2, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
dev.off()

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_boxplot-SEC_relativeIntensities.pdf'), width = 10, height = 6)
df_test %>%
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  ggplot(aes(x = response, y = intensity.relative, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~cell.type, scales = "free_y") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Comparison of cell type composition across response groups",
    subtitle = "Relative proportions of cell types based on EV proteins",
    x = "Response",
    y = "Relative proportion",
    fill = "Response",
    color = "Response"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = -0.5, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )
dev.off()

#### Beta regression
cell_types_EV <- colnames(EV_data[,c(3,4,5)])

results_beta_EV <- list()

for (cell in cell_types_EV) {
  
  y <- EV_data[[cell]] / 100   # percentá → proporcie
  
  # ošetrenie 0 a 1
  y[y == 0] <- 1e-6
  y[y == 1] <- 1 - 1e-6
  
  model_EV <- betareg(y ~ response, data = EV_data)
  
  # emmeans
  emm_EV <- emmeans(model_EV, ~ response)
  pairwise_EV <- pairs(emm_EV, adjust = "fdr")
  
  results_beta_EV[[cell]] <- list(
    model_EV = model_EV,
    emmeans_EV = emm_EV,
    pairs_EV = pairwise_EV
  )
}

beta_pvals_EV <- sapply(results_beta_EV, function(x) {
  anova(x$model_EV)$`Pr(>F)`[1]
})

############## Extended cell types (from flow cytometry analysis) #################

df_all_cells <- read_excel(here("data", "flowcytometry_data.xlsx"))

df_all_cells <- df_all_cells %>%
  rename(patient_id = `Patient ID`) %>%
  left_join(id_map, by = c("patient_id" = "old_id"))

# barplot for all patients
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_rawIntensities_all_cells.pdf'), width = 38, height = 15)
df_all_cells %>%
  pivot_longer(
    cols = c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
             "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells"),
    names_to = "cell.type",
    values_to = "proportion"
  ) %>%
  
  # prepočet percent 
  mutate(proportion = proportion / 100) %>%
  
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
    "Other cells" = "black"   # 👉 doplnená farba
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
    axis.text = element_text(size = 38), # čísla + sample_factor 
    axis.text.x = element_text(angle = 45, hjust = 1), # otočenie + veľkosť 
    legend.title = element_text(size = 40), # legenda názov 
    legend.text = element_text(size = 38),
    legend.position = "right"
  )
dev.off()

### better layout
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_heatmap-SEC_relativeIntensities_response_all_cells.pdf'), width = 30, height = 15)
# 📊 príprava dát
df_long <- df_all_cells %>%
  pivot_longer(
    cols = c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
             "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells"),
    names_to = "cell.type",
    values_to = "proportion"
  ) %>%
  mutate(proportion = proportion / 100) %>%
  mutate(patient_factor = factor(new_id, levels = sort(unique(new_id))))

# 🎨 spoločný theme
common_theme <- theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 60, face = "bold", hjust = 0),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 38),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 40),
    legend.text = element_text(size = 38),
    strip.text = element_text(size = 40, face = "bold")
  )

# 🔝 horný graf (RP + RT)
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
    "Other cells" = "black"   # 👉 doplnená farba
  )) +
  labs(
    x = "Patient",
    y = "Relative proportion",
    fill = "Cell Type"
  ) +
  common_theme

# 🔻 dolný graf (SE cez celú šírku)
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
    "Other cells" = "black"   # 👉 doplnená farba
  )) +
  labs(
    x = "Patient",
    y = "Relative proportion",
    fill = "Cell Type"
  ) +
  ggtitle("SE") +
  common_theme +
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5)
  )

# 🧩 spojenie (RP RT / SE SE)
final_plot <- p_top / p_bottom +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "right")

# ➕ titulok
final_plot +
  plot_annotation(
    title = "Relative cell type composition across patients by response",
    subtitle = "Based on FC data",
    theme = theme(
      plot.title = element_text(size = 55, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 40, face = "bold")
    )
  )
dev.off()

########################################################
################# boxploty + Wilcoxon ##################
################## Transformed data ####################
########################################################

pseudocount <- 1e-6

# 🔄 CLR transformácia
df_clr <- df_long %>%
  
  # každá vzorka = jedna kompozícia
  group_by(new_id) %>%   # <- zmeň na názov stĺpca so sample ID
  
  mutate(proportion = proportion + pseudocount) %>%
  
  # CLR
  mutate(
    clr_value = log(proportion / exp(mean(log(proportion))))
  ) %>%
  
  ungroup()

# 🔬 párové porovnania
comparisons <- list(
  c("RP", "RT"),
  c("RP", "SE"),
  c("RT", "SE")
)

pdf(
  here(
    'outputs',
    '5_RNAseq-comparison',
    '5_RNAseq-comparison_boxplot-SEC_CLR_all_cells.pdf'
  ),
  width = 15,
  height = 12
)

# 📊 boxplot CLR dát
df_clr %>%
  
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  
  ggplot(aes(x = response, y = clr_value, fill = response)) +
  
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  
  geom_jitter(
    aes(color = response),
    width = 0.1,
    alpha = 0.5,
    size = 1.5
  ) +
  
  facet_wrap(~cell.type, scales = "free_y", nrow = 2, ncol = 5) +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE,
    size=5
  ) +
  
  scale_fill_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  )) +
  
  scale_color_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  )) +
  
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

########################################################
################# boxploty + Wilcoxon ##################
################ Non-transformed data ##################
########################################################

df_test_cells <- df_all_cells %>%
  pivot_longer(
    cols = c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
             "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells"),
    names_to = "cell.type",
    values_to = "proportion"
  ) %>%
  
  # percentá → 0–1
  mutate(proportion = proportion / 100)

# 🔬 párové porovnania
comparisons <- list(
  c("RP", "RT"),
  c("RP", "SE"),
  c("RT", "SE")
)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_boxplot-SEC_percentages_all_cells.pdf'), width = 20, height = 23)
# 📦 boxplot
df_test_cells %>%
  mutate(response = factor(response, levels = c("SE", "RP", "RT"))) %>%
  ggplot(aes(x = response, y = proportion, fill = response)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = response), width = 0.1, alpha = 0.5, size = 1.5) +
  
  facet_wrap(~cell.type, scales = "free_y") +
  
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    p.adjust.method = "BH",
    hide.ns = TRUE
  ) +
  
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  
  labs(
    title = "Comparison of cell type composition across response groups",
    subtitle = "Relative proportions",
    x = "Response",
    y = "Relative proportion of cell types",
    fill = "Response",
    color = "Response"
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10),
    panel.spacing = unit(2, "lines"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

dev.off()


##### Beta regression #####
cell_types_all <- colnames(df_all_cells)
cell_types_all <- cell_types_all[-c(1,2,3)]

pvals <- c()

for (cell in cell_types_all) {
  y <- df_all_cells[[cell]] / 100
  
  # ošetrenie 0/1
  y[y == 0] <- 1e-6
  y[y == 1] <- 1 - 1e-6
  
  model <- betareg(y ~ df_all_cells$response)
  pvals[cell] <- summary(model)$coefficients$mean[2,4]
}

summary(model)

pvals_adj <- p.adjust(pvals, method = "fdr")

pvals_adj


emm <- emmeans(model, ~ response)
res <- pairs(emm, adjust = "fdr")
res

plot(emm)
emm_df <- as.data.frame(emm)

ggplot(emm_df, aes(x = response, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  ylim(0, NA) +
  labs(y = "Estimated proportion (model)", x = "Response group") 

###################################################
##################### rmcorr ######################
###################################################

#cells data
#df_cells_3 <- read_excel(here("data", "cells_percentages_3.xlsx"))
df_cells_3 <- df_all_cells %>%
  select(new_id, response, "Malignant cells", Fibroblasts, Macrophages)

df_cells_3 <- df_cells_3 %>%
  rename(
    malignant = `Malignant cells`,
    fibroblast = Fibroblasts,
    macrophage = Macrophages,
    patient_id = new_id
  )

df_cells_3_norm <- df_cells_3 %>%
  mutate(
    total = malignant + fibroblast + macrophage,
    
    malignant = malignant / total * 100,
    fibroblast = fibroblast / total * 100,
    macrophage = macrophage / total * 100
  ) %>%
  select(-total)


df_cells_3_long <- df_cells_3_norm %>%
  pivot_longer(
    cols = c(macrophage, fibroblast, malignant),
    names_to = "cell.type",
    values_to = "intensity.relative"
  ) %>%
  
  arrange(patient_id, response, cell.type)

df_cells_3_long <- df_cells_3_norm %>%
  pivot_longer(
    cols = c(macrophage, fibroblast, malignant),
    names_to = "cell.type",
    values_to = "FC_raw"
  ) %>%
  arrange(patient_id, response, cell.type)

# 🔄 CLR transformácia FC
df_cells_3_long <- df_cells_3_long %>%
  group_by(patient_id) %>%
  mutate(
    FC_raw = FC_raw + 1e-6,
    FC = log(FC_raw / exp(mean(log(FC_raw))))
  ) %>%
  ungroup()


# EVs data
df_test_percentages <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  mutate(
    patient_id = as.numeric(patient_id)
  ) %>%
  left_join(id_map, by = c("patient_id" = "old_id")) %>%
  group_by(new_id.x, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  
  group_by(new_id.x, response) %>%
  mutate(
    intensity.relative = mean.intensity / sum(mean.intensity),
    intensity.percent = intensity.relative * 100
  ) %>%
  
  ungroup() %>%
  
  arrange(new_id.x, response, cell.type)

df_test_percentages_long <- df_test_percentages[,c(1,2,3,6)]

df_test_percentages_long$patient_id <- as.numeric(df_test_percentages_long$patient_id)
df_test_percentages_long <- df_test_percentages_long %>%
  rename(patient_id = new_id.x)


## 
df_test_percentages <- SEC_long_intensity %>%
  separate(sample, into = c("patient_id", "response"), sep = "_") %>%
  mutate(patient_id = as.numeric(patient_id)) %>%
  left_join(id_map, by = c("patient_id" = "old_id")) %>%
  group_by(new_id.x, response, cell.type) %>%
  summarise(mean.intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  group_by(new_id.x, response) %>%
  mutate(
    intensity.relative = mean.intensity / sum(mean.intensity)
  ) %>%
  ungroup()

df_test_percentages_long <- df_test_percentages %>%
  select(new_id.x, response, cell.type, intensity.relative) %>%
  rename(patient_id = new_id.x)


# 🔄 CLR transformácia EV
df_test_percentages_long <- df_test_percentages_long %>%
  group_by(patient_id) %>%
  mutate(
    intensity.relative = intensity.relative + 1e-6,
    EV = log(intensity.relative / exp(mean(log(intensity.relative))))
  ) %>%
  ungroup()

# dataset left join for rmcorr
common_ids <- intersect(df_cells_3_long$patient_id,
                        df_test_percentages_long$patient_id)

df1 <- df_cells_3_long %>%
  filter(patient_id %in% common_ids)

df2 <- df_test_percentages_long %>%
  filter(patient_id %in% common_ids)

df_rmcorr <- left_join(
  df2,
  df1,
  by = c("patient_id", "cell.type")
)
df_rmcorr <- df_rmcorr %>%
  rename(FC = intensity.relative)
df_rmcorr <- df_rmcorr %>%
  rename(EV = intensity.percent)

## 
common_ids <- intersect(df_cells_3_long$patient_id,
                        df_test_percentages_long$patient_id)

df1 <- df_cells_3_long %>%
  filter(patient_id %in% common_ids)

df2 <- df_test_percentages_long %>%
  filter(patient_id %in% common_ids)

df_rmcorr <- left_join(
  df2,
  df1,
  by = c("patient_id", "cell.type")
)

# rmcorr
my.rmc.cells.EV <- rmcorr(participant = cell.type, measure1 = FC, measure2 = EV,
                           dataset = df_rmcorr)
df_rmcorr$fitted <- my.rmc.cells.EV$model$fitted.values

##
df_rmcorr %>%
  group_by(cell.type) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    n_obs = n()
  )

rmcorr_results <- df_rmcorr %>%
  group_by(cell.type) %>%
  group_map(~ rmcorr(
    participant = patient_id,
    measure1 = FC,
    measure2 = EV,
    dataset = .x
  ))

my.rmc.cells.EV <- rmcorr(
  participant = cell.type,
  measure1 = FC_raw,
  measure2 = intensity.relative,
  dataset = df_rmcorr
)

df_rmcorr$fitted <- my.rmc.cells.EV$model$fitted.values

pred_data <- df_rmcorr %>%
  group_by(cell.type) %>%
  group_map(~ {
    
    df <- filter(.x,
                 !is.na(FC_raw),
                 !is.na(intensity.relative))
    
    model <- lm(intensity.relative ~ FC_raw, data = df)
    
    x_seq <- seq(min(df$FC_raw), max(df$FC_raw), length.out = 100)
    
    data.frame(
      cell.type = unique(df$cell.type),
      FC_raw = x_seq,
      intensity.relative = predict(model,
                                   newdata = data.frame(FC_raw = x_seq))
    )
  }) %>%
  bind_rows()

# rmcorr plot
pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_rmcorr_FCvsEV.pdf'), width = 10, height = 8)
ggplot(df_rmcorr, aes(x = FC, y = EV, color = cell.type, group = cell.type)) + 
  geom_line(aes(x = FC, y = my.rmc.cells.EV$model$fitted.values), size = 1) + 
  geom_point() + 
  theme_prism() + 
  xlim(0,100) +
  ylim(0, 100) +
  labs(x = "FC (%)", y = "EV (%)", title = "Repeated measure correlation") +
  scale_color_manual(values = c(
    "malignant"="#F8766D",
    "fibroblast"="#00BA38",
    "macrophage"="#619CFF"
  ))+
  theme(
    plot.title = element_text(size = 25, hjust = 0.5), # názov grafu
    axis.title = element_text(size = 24),              # názvy osí
    axis.text = element_text(size = 22),               # hodnoty na osiach
    legend.title = element_text(size = 15),            # názov legendy
    legend.text = element_text(size = 13)              # položky legendy
  ) +
  annotate("text", x = 60, y = 14, label = paste("r", "=", round(my.rmc.cells.EV$r, digits = 3), sep = " "), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells.EV$CI[1], digits = 3), "-", round(my.rmc.cells.EV$CI[2], digits = 3)), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 2, label = paste("p = ", format(my.rmc.cells.EV$p, scientific = TRUE, digits = 3)), hjust = 0, size = 5) 
dev.off()

##
slopes_df <- df_rmcorr %>%
  group_by(cell.type) %>%
  summarise(
    slope = coef(lm(intensity.relative * 100 ~ FC_raw))[2],
    intercept = coef(lm(intensity.relative * 100 ~ FC_raw))[1]
  )

corr_df <- df_rmcorr %>%
  group_by(cell.type) %>%
  summarise(
    r = cor(FC_raw, intensity.relative * 100, method = "spearman", use = "complete.obs"),
    p = cor.test(FC_raw, intensity.relative * 100)$p.value
  ) %>%
  mutate(
    r_label = case_when(
      cell.type == "fibroblast" ~ paste0("italic(r)[f] == ", round(r, 2)),
      cell.type == "malignant"  ~ paste0("italic(r)[m] == ", round(r, 2)),
      cell.type == "macrophage" ~ paste0("italic(r)[mac] == ", round(r, 2))
    )
  )

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_rmcorr_FCvsEV_new.pdf'), width = 10, height = 8)
ggplot(df_rmcorr, aes(x = FC_raw, y = intensity.relative*100, color = cell.type)) +
  # Vykreslenie surových bodov
  geom_point(alpha = 0.7) + 
  
  # Pridanie samostatných regresných čiar pre každý cell.type
  # se = FALSE odstráni šedý interval spoľahlivosti okolo čiar (ak ho chceš, daj TRUE)
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 1) +
  geom_text(
    data = corr_df,
    aes(
      x = 80,
      y = c(95, 88, 81),
      label = r_label,
      color = cell.type
    ),
    parse = TRUE,
    size = 5,
    show.legend = FALSE
  )+
  
  theme_prism() +
  xlim(0, 100) +
  ylim(0, 100) +
  labs(x = "FC (%)", y = "EV (%)", title = "Repeated measure correlation") +
  
  scale_color_manual(values = c(
    "malignant" = "#F8766D",
    "fibroblast" = "#00BA38",
    "macrophage" = "#619CFF"
  )) +
  
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  ) +
  
  # Anotácie pre celkový rmcorr model (zostávajú podľa pôvodného výpočtu)
  annotate("text", x = 60, y = 14, label = paste("r =", round(my.rmc.cells.EV$r, digits = 3)), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 8, label = paste("95% CI:", round(my.rmc.cells.EV$CI[1], digits = 3), "-", round(my.rmc.cells.EV$CI[2], digits = 3)), hjust = 0, size = 5) +
  annotate("text", x = 60, y = 2, label = paste("p =", format(my.rmc.cells.EV$p, scientific = TRUE, digits = 3)), hjust = 0, size = 5)
dev.off()

# separátne lineárne modely + predikcie
pred_data <- df_rmcorr %>%
  dplyr::group_by(cell.type) %>%
  dplyr::group_map(~ {
    
    df <- dplyr::filter(.x,
                        !is.na(FC_raw),
                        !is.na(intensity.relative))
    
    if (nrow(df) < 2) return(NULL)
    
    model <- lm(intensity.relative ~ FC_raw, data = df)
    
    x_seq <- seq(min(df$FC_raw), max(df$FC_raw), length.out = 100)
    
    data.frame(
      cell.type = unique(df$cell.type),
      FC_raw = x_seq,
      intensity.relative = predict(
        model,
        newdata = data.frame(FC_raw = x_seq)
      )
    )
  }) %>%
  dplyr::bind_rows()


# GGPLOT
ggplot(df_rmcorr,
       aes(x = FC_raw,
           y = intensity.relative,
           color = cell.type)) +
  
  geom_point() +
  
  geom_line(
    data = pred_data,
    aes(x = FC_raw,
        y = intensity.relative,
        color = cell.type),
    size = 1
  ) +
  
  theme_prism() +
  
  xlim(0, 100) +
  ylim(0, 100) +
  
  labs(
    x = "FC (%)",
    y = "EV (%)",
    title = "Repeated measure correlation"
  ) +
  
  scale_color_manual(values = c(
    "malignant" = "#F8766D",
    "fibroblast" = "#00BA38",
    "macrophage" = "#619CFF"
  )) +
  
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  ) +
  
  annotate(
    "text", x = 60, y = 14,
    label = paste("r =", round(my.rmc.cells.EV$r, 3)),
    hjust = 0, size = 5
  ) +
  
  annotate(
    "text", x = 60, y = 8,
    label = paste(
      "95% CI:",
      round(my.rmc.cells.EV$CI[1], 3),
      "-",
      round(my.rmc.cells.EV$CI[2], 3)
    ),
    hjust = 0, size = 5
  ) +
  
  annotate(
    "text", x = 60, y = 2,
    label = paste(
      "p =",
      format(my.rmc.cells.EV$p, scientific = TRUE, digits = 3)
    ),
    hjust = 0, size = 5
  )

ggplot(df_rmcorr,
       aes(x = FC_raw, y = intensity.relative, color = cell.type)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ cell.type)

df_rmcorr %>%
  group_by(cell.type) %>%
  summarise(
    n = n(),
    n_na_fc = sum(is.na(FC_raw)),
    n_na_ev = sum(is.na(intensity.relative))
  )

## 
ggplot(df_rmcorr, aes(x = FC_raw, y = intensity.relative, color = cell.type)) +
  geom_point() +
  geom_line(aes(y = fitted, group = cell.type), size = 1) +
  theme_prism() +
  labs(
    x = "FC (CLR)",
    y = "EV (CLR)",
    title = "Repeated measures correlation (CLR-transformed data)"
  ) +
  scale_color_manual(values = c(
    "malignant" = "#F8766D",
    "fibroblast" = "#00BA38",
    "macrophage" = "#619CFF"
  )) +
  theme(
    plot.title = element_text(size = 25, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  ) +
  annotate("text", x = 1, y = 0.2,
           label = paste0("r = ", round(my.rmc.cells.EV$r, 3)),
           hjust = 0, size = 5) +
  annotate("text", x = 1, y = 0.1,
           label = paste0("95% CI: ",
                          round(my.rmc.cells.EV$CI[1], 3), " - ",
                          round(my.rmc.cells.EV$CI[2], 3)),
           hjust = 0, size = 5) +
  annotate("text", x = 1, y = 0,
           label = paste0("p = ", format(my.rmc.cells.EV$p, scientific = TRUE)),
           hjust = 0, size = 5)

# toto nefunguje - lmer potrebuje outcome (response) ako numeric
m1_EV <- lmer(response ~ EV + (1|patient_id), data = df_rmcorr)
m2_FC <- lmer(response ~ FC + (1|patient_id), data = df_rmcorr)

install.packages("ordinal")
library(ordinal)

df_rmcorr$response <- factor(df_rmcorr$response, 
                             levels = c("SE", "RP", "RT"),
                             ordered = TRUE)

clmm(response ~ EV + (1|patient_id), data = df_rmcorr)

###################################################################
############### clr + permanova (FC all cell types) ###############
###################################################################

cell_data <- df_all_cells[, c("Malignant cells", "Fibroblasts", "Macrophages", "Monocytes",
                    "B cells", "CD4 T cells", "CD8 T cells", "Neutrophils", "NK cells", "Other cells")]

# pseudocount
cell_data[cell_data == 0] <- 1e-6

# clr transformation
clr_data <- clr(cell_data)
dist_matrix <- dist(clr_data)

# permanova
adonis2(dist_matrix ~ response, data = df_all_cells)

# PCA
pca <- prcomp(clr_data)

pca_df <- data.frame(pca$x, response = df_all_cells$response)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_all_cells_raw.pdf'), width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = response)) +
  scale_color_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  ))+
  geom_point(size = 6) +
  theme(
    axis.title = element_text(size = 35),   # názvy osí (PC1, PC2)
    axis.text = element_text(size = 29),    # hodnoty na osiach
    legend.title = element_text(size = 31), # názov legendy
    legend.text = element_text(size = 29)   # položky v legende
  )
dev.off()

########################################################
########### clr + permanova (FC 3 cell types) ##########
########################################################

cell_data.3 <- df_cells_3_norm[,c("malignant", "fibroblast", "macrophage")]

# pseudocount
cell_data.3[cell_data.3 == 0] <- 1e-6

# clr transformation
clr_data.3 <- clr(cell_data.3)
dist_matrix.3 <- dist(clr_data.3)

# permanova
adonis2(dist_matrix.3 ~ response, data = df_cells_3_norm)

# PCA
pca.3 <- prcomp(clr_data.3)

pca_df.3 <- data.frame(pca.3$x, response = df_cells_3_norm$response)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_cells_3_raw.pdf'), width = 10, height = 8)
ggplot(pca_df.3, aes(x = PC1, y = PC2, color = response)) +
  geom_point(size = 6) +
  scale_color_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  ))+
  theme(
    axis.title = element_text(size = 35),   # názvy osí (PC1, PC2)
    axis.text = element_text(size = 29),    # hodnoty na osiach
    legend.title = element_text(size = 31), # názov legendy
    legend.text = element_text(size = 29)   # položky v legende
  )
dev.off()

####################################################
############## clr + permanova (EV) ################
####################################################

EV_data <- df_test_percentages_long %>%
  pivot_wider(
    names_from = cell.type,
    values_from = intensity.percent
  )

EV_data_cells <- EV_data[,c("malignant", "fibroblast", "macrophage")]
EV_data_cells[EV_data_cells == 0] <- 1e-6

clr_data_EV <- clr(EV_data_cells)
dist_matrix_EV <- dist(clr_data_EV)

adonis2(dist_matrix_EV ~ response, data = EV_data)

# PCA
pca_EV <- prcomp(clr_data_EV)
pca_df_EV <- data.frame(pca_EV$x, response = EV_data$response)

pdf(here('outputs', '5_RNAseq-comparison', '5_RNAseq-comparison_PCA_EVs_raw.pdf'), width = 10, height = 8)
ggplot(pca_df_EV, aes(x = PC1, y = PC2, color = response)) +
  geom_point(size = 6) +
  scale_color_manual(values = c(
    "RT" = "#D55E00",   # burnt orange
    "RP" = "#CC79A7",   # muted magenta
    "SE" = "#0072B2"    # deep blue
  ))+
  theme(
    axis.title = element_text(size = 35),   # názvy osí (PC1, PC2)
    axis.text = element_text(size = 29),    # hodnoty na osiach
    legend.title = element_text(size = 31), # názov legendy
    legend.text = element_text(size = 29)   # položky v legende
  )
dev.off()
