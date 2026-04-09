# =========================================================================
# Script Name: 04_Extract overlapping genes.R
# Description: Calculate the intersection between WGCNA module genes and 
#              Differentially Expressed Genes (DEGs), and generate a 
#              publication-ready customized Venn diagram.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

library(VennDiagram)
library(grid)

# Set working directory (Update this to your local repository path)
# setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\04_Extract overlapping genes")

# ====================== Data Loading ======================
# Read gene lists from text files (Assuming one gene per line, no header needed)
gene_set_A <- read.table("WGCNA.txt", header = FALSE, stringsAsFactors = FALSE)[,1]
gene_set_B <- read.table("DEGs.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

# ====================== Calculate Overlap ======================
# Calculate intersection and unique gene counts
overlap <- length(intersect(gene_set_A, gene_set_B))
only_A <- length(setdiff(gene_set_A, gene_set_B))
only_B <- length(setdiff(gene_set_B, gene_set_A))

# ====================== Plotting Venn Diagram ======================
pdf("Venn_red_blue_horizontal.pdf", width = 6, height = 6)

# Draw Venn diagram circles (disabling default number display for custom alignment)
venn.plot <- draw.pairwise.venn(
  area1 = 100,           # Equal circle sizes for aesthetic symmetry
  area2 = 100,
  cross.area = 50,       # Placeholder for intersection area to enforce overlap shape
  category = c("WGCNA", "DEGs"),
  lty = 1,
  lwd = 2,
  col = "black",
  fill = c("#D25E5E", "#4E8AC8"), # Red and Blue custom palette
  cat.col = "black",
  cat.cex = 1.4,
  cat.pos = c(-20, 20),
  cex = 0,               # ⭐ Disable default number rendering
  fontface = "bold",
  scaled = FALSE
)

# Optimize text positioning
# Manually place numbers using grid coordinates to ensure perfect centering
grid.text(label = only_A, x = 0.2, y = 0.48, gp = gpar(fontsize = 16, fontface = "bold"))
grid.text(label = only_B, x = 0.8, y = 0.48, gp = gpar(fontsize = 16, fontface = "bold"))
grid.text(label = overlap, x = 0.5, y = 0.48, gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()
