# =============================================================================
# Pseudotime Analysis using Monocle 2
# =============================================================================
# =============================================================================
# Pseudotime Analysis using Monocle 2
# Author: Hanzhang
# Description: Using monocle to ananlysis the pesudo time of postnatal liver 
#              development process and visualization
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggridges)

# =============================================================================
# Data Loading and Preprocessing
# =============================================================================

setwd("/path/to/your/working/directory")
load("phase.RData")
seu <- datH.comb

seu <- seu[, seu$seurat_clusters %in% c(0, 1, 2)]
seu <- RenameIdents(seu, `0` = "core", `1` = "late", `2` = "early")

DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5)

seu$cell_type <- Idents(seu)
table(seu$cell_type)

# =============================================================================
# CellDataSet Object Preparation
# =============================================================================

umi <- GetAssayData(seu, assay = "RNA", slot = "counts")

pData <- as.data.frame(seu@meta.data)
pData$celltype <- seu@active.ident

fData <- data.frame(
  gene_short_name = row.names(seu),
  row.names = row.names(seu)
)

common_genes <- intersect(row.names(umi), row.names(fData))

umi <- umi[common_genes, ]
fData <- fData[common_genes, ]

pd <- new('AnnotatedDataFrame', data = pData)
fd <- new('AnnotatedDataFrame', data = fData)

# =============================================================================
# Monocle 2 CellDataSet Creation
# =============================================================================

cds_monocle2 <- newCellDataSet(
  umi,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = negbinomial.size()
)

cds_monocle2 <- estimateSizeFactors(cds_monocle2)
cds_monocle2 <- estimateDispersions(cds_monocle2)

# =============================================================================
# Gene Detection and Filtering
# =============================================================================

cds_monocle2 <- detectGenes(cds_monocle2, min_expr = 0.1)

expressed_genes <- row.names(
  subset(Biobase::featureData(cds_monocle2)@data, 
         num_cells_expressed >= 10)
)

express_genes <- VariableFeatures(seu)
cds_monocle2 <- setOrderingFilter(cds_monocle2, express_genes)
plot_ordering_genes(cds_monocle2)

# =============================================================================
# Dimensionality Reduction and Trajectory Construction
# =============================================================================

cds_monocle2 <- reduceDimension(cds_monocle2, max_components = 2, method = 'DDRTree')

cds_monocle2 <- orderCells(cds_monocle2)
plot_cell_trajectory(cds_monocle2)

cds_monocle2 <- orderCells(cds_monocle2, root_state = 3)
plot_cell_trajectory(cds_monocle2)
plot_cell_trajectory(cds_monocle2, color_by = "Pseudotime")
plot_cell_trajectory(cds_monocle2, color_by = "seurat_clusters")

plot_cell_trajectory(cds_monocle2, color_by = "seurat_clusters") +
  facet_wrap("~seurat_clusters", nrow = 3)

plot_cell_trajectory(cds_monocle2, color_by = "cell_type") +
  facet_wrap("~cell_type", nrow = 3)

plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 3)

plot_complex_cell_trajectory(cds_monocle2, x = 1, y = 2, color_by = "seurat_clusters")

plot_cell_trajectory(cds_monocle2, color_by = "day") + 
  facet_wrap("~day", nrow = 3)

# =============================================================================
# Enhanced Trajectory Visualization
# =============================================================================

plot_trajectory <- plot_cell_trajectory(cds_monocle2, color_by = "day") +  
  facet_wrap(~ day, nrow = 2) +
  labs(x = "Component 1", y = "Component 2", color = "day") +
  scale_color_manual(values = c("P07" = "#440154FF",
                                "P14" = "#31688EFF",
                                "P21" = "#35B779FF",
                                "P28" = "#FDE725FF"))

print(plot_trajectory)

# =============================================================================
# Pseudotime Density Analysis
# =============================================================================

pseudotime_df <- as.data.frame(pData(cds_monocle2)[, c("Pseudotime", "day")])
pseudotime_df$cell <- rownames(pData(cds_monocle2))

head(pseudotime_df)

pseudotime_df <- pseudotime_df[order(pseudotime_df$Pseudotime), ]

p <- ggplot(pseudotime_df, aes(x = Pseudotime, fill = day, color = day)) +    
  geom_density(alpha = 0.5, size = 1.2) +  
  scale_fill_manual(values = c("P07" = "#440154FF",
                               "P14" = "#31688EFF",
                               "P21" = "#35B779FF",
                               "P28" = "#FDE725FF")) +
  scale_color_manual(values = c("P07" = "#440154FF",
                                "P14" = "#31688EFF",
                                "P21" = "#35B779FF",
                                "P28" = "#FDE725FF")) +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line = element_line(color = 'black'),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = 'black', size = 10)
  ) + 
  labs(
    title = "Pseudotime Density Plot by Day",
    x = "Pseudotime",
    y = "Density"
  )

print(p)

# =============================================================================
# Ridge Plot Visualization
# =============================================================================

p_ridge <- ggplot(pseudotime_df, aes(x = Pseudotime, y = day, fill = day)) +  
  geom_density_ridges(scale = 1) +
  scale_y_discrete(position = 'left') +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = 'black'),      
    axis.title.x = element_text(colour = 'black', size = 12),
    axis.title.y = element_text(colour = 'black', size = 12),
    axis.text = element_text(colour = 'black', size = 8)
  ) +
  scale_x_continuous(position = 'bottom') +
  scale_fill_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")) +
  labs(x = "Pseudotime", y = "Day")

print(p_ridge)
