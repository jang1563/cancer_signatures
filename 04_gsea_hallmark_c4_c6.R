# =============================================================================
# GSEA Analysis - Hallmark, C4, and C6 Gene Sets (I4 Study)
# =============================================================================
#
# Purpose: Visualize gene set enrichment analysis (fGSEA) results for PBMC
#          and immune cell subsets comparing post-flight to pre-flight conditions
#
# Author: [Your Name]
# Date: April 2025
# Repository: https://github.com/[username]/cancer_signatures
#
# Input:
#   - Pre-computed fGSEA results (RDS files) for each cell type
#   - MSigDB gene sets: Hallmark, C4 (computational cancer), C6 (oncogenic)
#
# Output:
#   - Dot plots showing NES and significance across cell types and timepoints
#   - Exported CSV tables for supplementary materials
#
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(fgsea)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggthemes)
library(ggtext)
library(stringr)
library(forcats)
library(data.table)

# -----------------------------------------------------------------------------
# 2. Define Cell Types
# -----------------------------------------------------------------------------
cell_list <- c('pbmc', "CD4_T", "CD8_T", "other_T", "B", "NK",
               "CD14_Mono", "CD16_Mono", "DC", 'other')

# -----------------------------------------------------------------------------
# 3. Load Pre-computed fGSEA Results
# -----------------------------------------------------------------------------
for (i in cell_list) {
  # Hallmark gene sets
  a4 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.hallmark.immediately.rds"
  ))
  a5 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.hallmark.longterm.rds"
  ))

  # C4 gene sets (computational cancer)
  b4 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.C4.immediately.rds"
  ))
  b5 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.C4.longterm.rds"
  ))

  # C6 gene sets (oncogenic signatures)
  c4 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.C6.immediately.rds"
  ))
  c5 <- readRDS(paste0(
    "../../Differential_expression/with_filter/Enrichment_analysis/",
    "2022_10_27_fgsea_Integration_plasma_exosome_skin/2022_10_30_PBMCs/2023_12_01/intermediate/",
    i, ".fgsea.C6.longterm.rds"
  ))

  # Assign to cell-type specific variables
  assign(paste0(i, ".fgsea.hallmark.immediately"), a4)
  assign(paste0(i, ".fgsea.hallmark.longterm"), a5)
  assign(paste0(i, ".fgsea.C4.immediately"), b4)
  assign(paste0(i, ".fgsea.C4.longterm"), b5)
  assign(paste0(i, ".fgsea.C6.immediately"), c4)
  assign(paste0(i, ".fgsea.C6.longterm"), c5)
}

# -----------------------------------------------------------------------------
# 4. Process and Combine Results
# -----------------------------------------------------------------------------
for (i in cell_list) {
  # Get data for each cell type
  a4 <- get(paste0(i, ".fgsea.hallmark.immediately"))
  a5 <- get(paste0(i, ".fgsea.hallmark.longterm"))
  b4 <- get(paste0(i, ".fgsea.C4.immediately"))
  b5 <- get(paste0(i, ".fgsea.C4.longterm"))
  c4 <- get(paste0(i, ".fgsea.C6.immediately"))
  c5 <- get(paste0(i, ".fgsea.C6.longterm"))

  # Add timepoint labels
  a4$timepoint <- "R+1"
  a5$timepoint <- "R+45&R+82"
  b4$timepoint <- "R+1"
  b5$timepoint <- "R+45&R+82"
  c4$timepoint <- "R+1"
  c5$timepoint <- "R+45&R+82"

  # Add category labels
  a4$info <- "Hallmark"
  a5$info <- "Hallmark"
  b4$info <- "C4"
  b5$info <- "C4"
  c4$info <- "C6"
  c5$info <- "C6"

  # Add cell type labels
  a4$celltype <- i
  a5$celltype <- i
  b4$celltype <- i
  b5$celltype <- i
  c4$celltype <- i
  c5$celltype <- i

  # Combine Hallmark
  hallmark_combined <- rbind(a4, a5)
  assign(paste0(i, ".hallmark.2"), hallmark_combined)

  # Combine C4 (separate up/down for visualization)
  c4_up <- b4 %>% filter(NES > 0)
  c4_down <- b4 %>% filter(NES < 0)
  c4_long_up <- b5 %>% filter(NES > 0)
  c4_long_down <- b5 %>% filter(NES < 0)
  assign(paste0(i, ".C4.2.up"), c4_up)
  assign(paste0(i, ".C4.2.down"), c4_down)
  assign(paste0(i, ".C4.long.up"), c4_long_up)
  assign(paste0(i, ".C4.long.down"), c4_long_down)

  # Combine C6
  c6_up <- c4 %>% filter(NES > 0)
  c6_down <- c4 %>% filter(NES < 0)
  c6_long_up <- c5 %>% filter(NES > 0)
  c6_long_down <- c5 %>% filter(NES < 0)
  assign(paste0(i, ".C6.2.up"), c6_up)
  assign(paste0(i, ".C6.2.down"), c6_down)
  assign(paste0(i, ".C6.long.up"), c6_long_up)
  assign(paste0(i, ".C6.long.down"), c6_long_down)
}

# -----------------------------------------------------------------------------
# 5. Combine All Cell Types
# -----------------------------------------------------------------------------
hallmark.2 <- rbind(
  pbmc.hallmark.2, CD4_T.hallmark.2, CD8_T.hallmark.2, other_T.hallmark.2,
  B.hallmark.2, NK.hallmark.2, CD14_Mono.hallmark.2, CD16_Mono.hallmark.2,
  DC.hallmark.2, other.hallmark.2
)

C4.2 <- rbind(
  pbmc.C4.2.up, pbmc.C4.2.down, CD4_T.C4.2.up, CD4_T.C4.2.down,
  CD8_T.C4.2.up, CD8_T.C4.2.down, other_T.C4.2.up, other_T.C4.2.down,
  B.C4.2.up, B.C4.2.down, NK.C4.2.up, NK.C4.2.down,
  CD14_Mono.C4.2.up, CD14_Mono.C4.2.down, CD16_Mono.C4.2.up, CD16_Mono.C4.2.down,
  DC.C4.2.up, DC.C4.2.down, other.C4.2.up, other.C4.2.down,
  pbmc.C4.long.up, pbmc.C4.long.down, CD4_T.C4.long.up, CD4_T.C4.long.down,
  CD8_T.C4.long.up, CD8_T.C4.long.down, other_T.C4.long.up, other_T.C4.long.down,
  B.C4.long.up, B.C4.long.down, NK.C4.long.up, NK.C4.long.down,
  CD14_Mono.C4.long.up, CD14_Mono.C4.long.down, CD16_Mono.C4.long.up, CD16_Mono.C4.long.down,
  DC.C4.long.up, DC.C4.long.down, other.C4.long.up, other.C4.long.down
)

C6.2 <- rbind(
  pbmc.C6.2.up, pbmc.C6.2.down, CD4_T.C6.2.up, CD4_T.C6.2.down,
  CD8_T.C6.2.up, CD8_T.C6.2.down, other_T.C6.2.up, other_T.C6.2.down,
  B.C6.2.up, B.C6.2.down, NK.C6.2.up, NK.C6.2.down,
  CD14_Mono.C6.2.up, CD14_Mono.C6.2.down, CD16_Mono.C6.2.up, CD16_Mono.C6.2.down,
  DC.C6.2.up, DC.C6.2.down, other.C6.2.up, other.C6.2.down,
  pbmc.C6.long.up, pbmc.C6.long.down, CD4_T.C6.long.up, CD4_T.C6.long.down,
  CD8_T.C6.long.up, CD8_T.C6.long.down, other_T.C6.long.up, other_T.C6.long.down,
  B.C6.long.up, B.C6.long.down, NK.C6.long.up, NK.C6.long.down,
  CD14_Mono.C6.long.up, CD14_Mono.C6.long.down, CD16_Mono.C6.long.up, CD16_Mono.C6.long.down,
  DC.C6.long.up, DC.C6.long.down, other.C6.long.up, other.C6.long.down
)

# Combine all datasets
pbmc.sub.2 <- rbind(hallmark.2, C4.2, C6.2)

# -----------------------------------------------------------------------------
# 6. Clean Cell Type Labels
# -----------------------------------------------------------------------------
pbmc.sub.2$celltype <- ifelse(pbmc.sub.2$celltype == "pbmc", "PBMC",
                       ifelse(pbmc.sub.2$celltype == "CD4_T", "CD4 T",
                       ifelse(pbmc.sub.2$celltype == "CD8_T", "CD8 T",
                       ifelse(pbmc.sub.2$celltype == "other_T", "other T",
                       ifelse(pbmc.sub.2$celltype == "B", "B",
                       ifelse(pbmc.sub.2$celltype == "NK", "NK",
                       ifelse(pbmc.sub.2$celltype == "CD14_Mono", "CD14 Mono",
                       ifelse(pbmc.sub.2$celltype == "CD16_Mono", "CD16 Mono",
                       ifelse(pbmc.sub.2$celltype == "DC", "DC",
                       ifelse(pbmc.sub.2$celltype == "other", "other", pbmc.sub.2$celltype))))))))))

# Set factor levels for consistent ordering
pbmc.sub.2$celltype <- factor(pbmc.sub.2$celltype,
                              levels = c('PBMC', "CD4 T", "CD8 T", "other T", "B",
                                        "NK", "CD14 Mono", "CD16 Mono", "DC", "other"))

# -----------------------------------------------------------------------------
# 7. Create Output Directories
# -----------------------------------------------------------------------------
dir.create('intermediate/2025_01_17', recursive = TRUE, showWarnings = FALSE)
dir.create('intermediate/2025_04_24/', recursive = TRUE, showWarnings = FALSE)
dir.create('output/2025_01_17/', recursive = TRUE, showWarnings = FALSE)

# Save intermediate RDS
saveRDS(pbmc.sub.2, "intermediate/2025_01_17/pbmc.Hallmark.C4.C6.sub.2.rds")

# -----------------------------------------------------------------------------
# 8. Helper Function for Integer Axis Breaks
# -----------------------------------------------------------------------------
int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}

# -----------------------------------------------------------------------------
# 9. Define Plot Labels
# -----------------------------------------------------------------------------
Category.labs <- c("Hallmark", "C4", "C6")
names(Category.labs) <- c("Hallmark", "C4", "C6")

sample.labs <- c("PBMC", "CD4 T", "CD8 T", "other T", "B",
                 "NK", "CD14 Mono", "CD16 Mono", "DC", "other")
names(sample.labs) <- c("PBMC", "CD4 T", "CD8 T", "other T", "B",
                        "NK", "CD14 Mono", "CD16 Mono", "DC", "other")

# -----------------------------------------------------------------------------
# 10. Prepare Hallmark Data and Export
# -----------------------------------------------------------------------------
pbmc.Hallmark <- pbmc.sub.2 %>%
  dplyr::filter(info == 'Hallmark') %>%
  dplyr::filter(padj < 0.05)

# Export for supplementary table
df_export <- as.data.frame(pbmc.Hallmark)
list_cols <- vapply(df_export, is.list, logical(1))
df_export[, list_cols] <- lapply(df_export[, list_cols, drop = FALSE],
                                  function(x) sapply(x, paste, collapse = ", "))
write.csv(df_export, "intermediate/2025_04_24/Hallmark.I4.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Plot Hallmark Gene Sets
# -----------------------------------------------------------------------------
options(repr.plot.width = 24, repr.plot.height = 10)

p_hallmark <- ggplot(pbmc.Hallmark,
                     aes(x = timepoint,
                         y = fct_reorder(pathway, NES),
                         size = abs(-log10(padj)),
                         fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~celltype, nrow = 1, labeller = labeller(celltype = sample.labs)) +
  labs(
    title = "Hallmark Gene Sets (padj < 0.05)",
    x = "Timepoint",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave("output/2025_01_17/I4.Hallmark_dotplot.pdf", p_hallmark, width = 24, height = 10)
ggsave("output/2025_01_17/I4.Hallmark_dotplot.tiff", p_hallmark,
       width = 24, height = 10, dpi = 300, compression = 'lzw')

# -----------------------------------------------------------------------------
# 12. Prepare and Plot C4 Gene Sets
# -----------------------------------------------------------------------------
pbmc.C4 <- pbmc.sub.2 %>%
  filter(info == 'C4') %>%
  filter(padj < 0.001)

# Export for supplementary table
df_export <- as.data.frame(pbmc.sub.2 %>% filter(info == 'C4'))
list_cols <- vapply(df_export, is.list, logical(1))
df_export[, list_cols] <- lapply(df_export[, list_cols, drop = FALSE],
                                  function(x) sapply(x, paste, collapse = ", "))
write.csv(df_export, "intermediate/2025_04_24/C4.I4.csv", row.names = FALSE)

options(repr.plot.width = 24, repr.plot.height = 40)

p_c4 <- ggplot(pbmc.C4,
               aes(x = timepoint,
                   y = fct_reorder(pathway, NES),
                   size = abs(-log10(padj)),
                   fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~celltype, nrow = 1, labeller = labeller(celltype = sample.labs)) +
  labs(
    title = "C4 Gene Sets - Computational Cancer (padj < 0.001)",
    x = "Timepoint",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave("output/2025_01_17/I4.C4_dotplot.pdf", p_c4, width = 24, height = 40)
ggsave("output/2025_01_17/I4.C4_dotplot.tiff", p_c4,
       width = 24, height = 40, dpi = 300, compression = 'lzw')

# -----------------------------------------------------------------------------
# 13. Prepare and Plot C6 Gene Sets
# -----------------------------------------------------------------------------
pbmc.C6 <- pbmc.sub.2 %>%
  filter(info == 'C6') %>%
  filter(padj < 0.05)

# Export for supplementary table
df_export <- as.data.frame(pbmc.sub.2 %>% filter(info == 'C6'))
list_cols <- vapply(df_export, is.list, logical(1))
df_export[, list_cols] <- lapply(df_export[, list_cols, drop = FALSE],
                                  function(x) sapply(x, paste, collapse = ", "))
write.csv(df_export, "intermediate/2025_04_24/C6.I4.csv", row.names = FALSE)

options(repr.plot.width = 24, repr.plot.height = 20)

p_c6 <- ggplot(pbmc.C6,
               aes(x = timepoint,
                   y = fct_reorder(pathway, NES),
                   size = abs(-log10(padj)),
                   fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~celltype, nrow = 1, labeller = labeller(celltype = sample.labs)) +
  labs(
    title = "C6 Gene Sets - Oncogenic Signatures (padj < 0.05)",
    x = "Timepoint",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave("output/2025_01_17/I4.C6_dotplot.pdf", p_c6, width = 24, height = 20)
ggsave("output/2025_01_17/I4.C6_dotplot.tiff", p_c6,
       width = 24, height = 20, dpi = 300, compression = 'lzw')

# -----------------------------------------------------------------------------
# 14. Session Info
# -----------------------------------------------------------------------------
sessionInfo()
