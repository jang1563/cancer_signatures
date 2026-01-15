# =============================================================================
# GSEA Analysis - NASA Twins Study: Hallmark, C4, and C6 Gene Sets
# =============================================================================
#
# Purpose: Visualize gene set enrichment analysis (fGSEA) results comparing
#          in-flight and post-flight conditions to pre-flight baseline
#          using NASA Twins Study data
#
# Author: [Your Name]
# Date: January 2025
# Repository: https://github.com/[username]/cancer_signatures
#
# Input:
#   - NASA Twins Study fGSEA output (CSV)
#   - MSigDB gene sets: Hallmark, C4 (computational cancer), C6 (oncogenic)
#
# Output:
#   - Dot plots showing NES and significance across cell types and conditions
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
library(RColorBrewer)
library(data.table)
library(msigdbr)

# -----------------------------------------------------------------------------
# 2. Load NASA Twins GSEA Data
# -----------------------------------------------------------------------------
df <- read.csv("../../NASA_twinstudy/NASA_Twins_fGSEA_output/Twins.GSEA.csv")

# Create combined celltype identifiers
df$celltype <- paste0(df$CellType, ";", df$ExperimentType)
df$Celltype <- paste0(df$CellType, "_", df$ExperimentType)

# -----------------------------------------------------------------------------
# 3. Clean Cell Type Labels
# -----------------------------------------------------------------------------
df$Celltype <- ifelse(df$Celltype == 'CD4_Multivariate (PolyA+ and Ribodepleted together)', 'CD4_Multivariate',
               ifelse(df$Celltype == 'CD4_PolyA+', 'CD4_PolyA',
               ifelse(df$Celltype == 'CD4_Ribodepleted', 'CD4_Ribodepleted',
               ifelse(df$Celltype == 'CD8_Multivariate (PolyA+ and Ribodepleted together)', 'CD8_Multivariate',
               ifelse(df$Celltype == 'CD8_PolyA+', 'CD8_PolyA',
               ifelse(df$Celltype == 'CD8_Ribodepleted', 'CD8_Ribodepleted',
               ifelse(df$Celltype == 'LD_Multivariate (PolyA+ and Ribodepleted together)', 'LD_Multivariate',
               ifelse(df$Celltype == 'LD_PolyA+', 'LD_PolyA',
               ifelse(df$Celltype == 'LD_Ribodepleted', 'LD_Ribodepleted',
               ifelse(df$Celltype == 'LN_Multivariate (PolyA+ and Ribodepleted together)', 'LN_Multivariate',
               ifelse(df$Celltype == 'LN_PolyA+', 'LN_PolyA',
               ifelse(df$Celltype == 'LN_Ribodepleted', 'LN_Ribodepleted', df$Celltype))))))))))))

# Standardize column name
colnames(df)[7] <- 'padj'

# Filter for relevant comparisons
df <- df %>% filter(Coefficient %in% c('In-flight vs Pre-flight', 'Post-flight vs Pre-flight'))
colnames(df)[2] <- 'timepoint'

# Convert to data.table for efficient filtering
dt <- as.data.table(df)

# -----------------------------------------------------------------------------
# 4. Load MSigDB Gene Sets
# -----------------------------------------------------------------------------
m_df.hu <- msigdbr(species = "Homo sapiens")

# Filter by category
m_df.hu.Hallmark <- m_df.hu %>% filter(gs_cat %in% c("H"))
m_df.hu.C4 <- m_df.hu %>% filter(gs_cat %in% c("C4"))
m_df.hu.C6 <- m_df.hu %>% filter(gs_cat %in% c("C6"))

# Create gene set lists
m_list.hu.Hallmark <- m_df.hu.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C4 <- m_df.hu.C4 %>% split(x = .$gene_symbol, f = .$gs_name)
m_list.hu.C6 <- m_df.hu.C6 %>% split(x = .$gene_symbol, f = .$gs_name)

# Get gene set names
Hallmark <- names(m_list.hu.Hallmark)
C4 <- names(m_list.hu.C4)
C6 <- names(m_list.hu.C6)

# -----------------------------------------------------------------------------
# 5. Filter Data by Gene Set Category
# -----------------------------------------------------------------------------
Hallmark.dt <- dt[grep(paste(Hallmark, collapse = "|"), Pathway, ignore.case = TRUE)]
C4.dt <- dt[grep(paste(C4, collapse = "|"), Pathway, ignore.case = TRUE)]
C6.dt <- dt[grep(paste(C6, collapse = "|"), Pathway, ignore.case = TRUE)]

# -----------------------------------------------------------------------------
# 6. Create Output Directories
# -----------------------------------------------------------------------------
dir.create('output/2025_01_17/NASA/', recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 7. Define Plot Labels
# -----------------------------------------------------------------------------
sample.labs <- c(
  'CD4;\nMultivariate \n(PolyA+ and \nRibodepleted \ntogether)',
  'CD4;\nPolyA+',
  'CD4;\nRibodepleted',
  'CD8;\nMultivariate \n(PolyA+ and \nRibodepleted \ntogether)',
  'CD8;\nPolyA+',
  'CD8;\nRibodepleted',
  'LD;\nMultivariate \n(PolyA+ and \nRibodepleted \ntogether)',
  'LD;\nPolyA+',
  'LD;\nRibodepleted',
  'LN;\nMultivariate \n(PolyA+ and \nRibodepleted \ntogether)',
  'LN;\nPolyA+',
  'LN;\nRibodepleted'
)

names(sample.labs) <- c(
  'CD4_Multivariate', 'CD4_PolyA', 'CD4_Ribodepleted',
  'CD8_Multivariate', 'CD8_PolyA', 'CD8_Ribodepleted',
  'LD_Multivariate', 'LD_PolyA', 'LD_Ribodepleted',
  'LN_Multivariate', 'LN_PolyA', 'LN_Ribodepleted'
)

# -----------------------------------------------------------------------------
# 8. Helper Function for Integer Axis Breaks
# -----------------------------------------------------------------------------
int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}

# -----------------------------------------------------------------------------
# 9. Plot Hallmark Gene Sets
# -----------------------------------------------------------------------------
options(repr.plot.width = 35, repr.plot.height = 20)

p_hallmark <- ggplot(Hallmark.dt,
                     aes(x = timepoint,
                         y = fct_reorder(Pathway, NES),
                         size = abs(-log10(padj)),
                         fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~Celltype, nrow = 1, labeller = labeller(Celltype = sample.labs)) +
  labs(
    title = "Hallmark Gene Sets - NASA Twins Study",
    x = "Comparison",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave(
  "output/2025_01_17/NASA/Hallmark_dotplot.pdf",
  p_hallmark,
  width = 35,
  height = 20
)

ggsave(
  "output/2025_01_17/NASA/Hallmark_dotplot.tiff",
  p_hallmark,
  width = 35,
  height = 20,
  dpi = 300,
  compression = 'lzw'
)

# -----------------------------------------------------------------------------
# 10. Plot C4 Gene Sets (Computational Cancer Gene Sets)
# -----------------------------------------------------------------------------
# Filter for highly significant results
C4.dt.filter <- C4.dt %>% filter(padj < 0.001)

options(repr.plot.width = 35, repr.plot.height = 45)

p_c4 <- ggplot(C4.dt.filter,
               aes(x = timepoint,
                   y = fct_reorder(Pathway, NES),
                   size = abs(-log10(padj)),
                   fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~Celltype, nrow = 1, labeller = labeller(Celltype = sample.labs)) +
  labs(
    title = "C4 Gene Sets (Computational Cancer) - NASA Twins Study",
    subtitle = "Filtered for padj < 0.001",
    x = "Comparison",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave(
  "output/2025_01_17/NASA/C4_dotplot.pdf",
  p_c4,
  width = 35,
  height = 45
)

ggsave(
  "output/2025_01_17/NASA/C4_dotplot.tiff",
  p_c4,
  width = 35,
  height = 45,
  dpi = 300,
  compression = 'lzw'
)

# -----------------------------------------------------------------------------
# 11. Plot C6 Gene Sets (Oncogenic Signatures)
# -----------------------------------------------------------------------------
options(repr.plot.width = 35, repr.plot.height = 35)

p_c6 <- ggplot(C6.dt,
               aes(x = timepoint,
                   y = fct_reorder(Pathway, NES),
                   size = abs(-log10(padj)),
                   fill = NES)) +
  geom_point(shape = 21, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_fill_gradientn(
    colors = c("navyblue", "white", "firebrick3"),
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  scale_size_continuous(range = c(2, 8)) +
  facet_wrap(~Celltype, nrow = 1, labeller = labeller(Celltype = sample.labs)) +
  labs(
    title = "C6 Gene Sets (Oncogenic Signatures) - NASA Twins Study",
    x = "Comparison",
    y = "Pathway",
    fill = "NES",
    size = "-Log10(padj)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

ggsave(
  "output/2025_01_17/NASA/C6_dotplot.pdf",
  p_c6,
  width = 35,
  height = 35
)

ggsave(
  "output/2025_01_17/NASA/C6_dotplot.tiff",
  p_c6,
  width = 35,
  height = 35,
  dpi = 300,
  compression = 'lzw'
)

# -----------------------------------------------------------------------------
# 12. Session Info
# -----------------------------------------------------------------------------
sessionInfo()
