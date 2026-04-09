# =========================================================================
# Script Name: 05_Enrichment Analysis.R
# Description: Perform Gene Ontology (GO) and KEGG pathway enrichment 
#              analysis on intersecting genes, and generate customized, 
#              publication-ready bar plots.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

rm(list = ls())

# Set working directory (Update this to your local repository path)
setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\05_Enrichment Analysis") 

# ------------------------------------------------------------------------------
# 1. Package Loading
# ------------------------------------------------------------------------------
packages <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "dplyr", "stringr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Data Preparation (Modified: Strictly read existing file)
# ------------------------------------------------------------------------------
# Check if the input file exists, stop execution if missing
if(!file.exists("Overlap_Genes.txt")) {
  stop("Error: 'Overlap_Genes.txt' not found in the working directory. Please ensure the file exists.")
}

# Read intersecting genes
genes <- read.table("Overlap_Genes.txt", header = F)
gene_list_symbol <- unique(str_trim(genes$V1))

# Convert Gene Symbols to Entrez IDs for clusterProfiler
entrez <- bitr(gene_list_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list <- entrez$ENTREZID

# ------------------------------------------------------------------------------
# 3. Define Global Color Palette (Red-Blue-Green)
# ------------------------------------------------------------------------------
uni_colors <- c(
  Red   = "#D25E5E",  
  Blue  = "#4E8AC8",  
  Green = "#4DAF4A"   
)

# ==============================================================================
# Part A: GO Enrichment Plot
# ==============================================================================
message(">>> Generating GO Enrichment Plot...")
go_result <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "ALL",
                      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

if(!is.null(go_result)) {
  # Extract top 5 terms per ontology category (BP, CC, MF)
  go_df <- as.data.frame(go_result) %>%
    group_by(ONTOLOGY) %>%
    arrange(pvalue) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    arrange(factor(ONTOLOGY, levels = c("MF", "CC", "BP")), Count)
  
  go_df$Description <- factor(go_df$Description, levels = go_df$Description)
  
  # Assign colors to ontologies
  cols_GO <- c("BP" = uni_colors[["Red"]], "CC" = uni_colors[["Blue"]], "MF" = uni_colors[["Green"]])
  
  pA <- ggplot(go_df, aes(x = Count, y = Description, fill = ONTOLOGY)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = Count), hjust = -0.3, size = 3) +
    scale_fill_manual(values = cols_GO) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "GO Enrichment Analysis", x = "Gene Count", y = NULL) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      
      # [Customization 1] Left-align title to the panel border
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.title.position = "panel", 
      
      legend.position = c(0.85, 0.2), 
      legend.background = element_rect(fill = "transparent")
    )
  
  pdf("Figure_A_GO.pdf", width = 8, height = 7)
  print(pA)
  dev.off()
}

# ==============================================================================
# Part B: KEGG Enrichment Plot (Customized layout and categorization)
# ==============================================================================
message(">>> Generating KEGG Enrichment Plot...")
kegg_result <- enrichKEGG(gene = gene_list, organism = 'hsa', pvalueCutoff = 1)

if(!is.null(kegg_result)) {
  kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Extract top 16 pathways
  kegg_df <- as.data.frame(kegg_result) %>%
    arrange(pvalue) %>%
    slice_head(n = 16) 
  
  # 1. Define custom functional categories based on keywords
  kegg_df$Category <- case_when(
    grepl("TNF|NF-kappa|Toll|IL-|Th17|Chemokine|Cytokine|NOD|Complement|Fc epsilon", 
          kegg_df$Description, ignore.case = T) ~ "Immune & Inflammation",
    grepl("bacterial|infection|Salmonella|Legionellosis|Staphylococcus|Pertussis|Chagas|Hepatitis|Influenza|Malaria", 
          kegg_df$Description, ignore.case = T) ~ "Pathogen Infection",
    grepl("MAPK|PI3K|Ras|Rap1|cAMP|cGMP|Signaling|FoxO|HIF-1|JAK-STAT", 
          kegg_df$Description, ignore.case = T) ~ "Signal Transduction",
    grepl("Apoptosis|Necroptosis|Ferroptosis|Autophagy|Metabolism|Glycolysis|Fluid shear|Lipid|Atherosclerosis", 
          kegg_df$Description, ignore.case = T) ~ "Metabolism & Cellular Process",
    TRUE ~ "Related Pathways"
  )
  
  # 2. Set factor levels for rendering order (Bottom to Top in ggplot)
  # Goal: Top (Red) -> Middle (Blue) -> Bottom (Green)
  sort_levels <- c(
    "Pathogen Infection", "Related Pathways",               # Bottom (Green)
    "Signal Transduction", "Metabolism & Cellular Process", # Middle (Blue)
    "Immune & Inflammation"                                 # Top (Red)
  )
  kegg_df$Category <- factor(kegg_df$Category, levels = sort_levels)
  
  # 3. Sort dataframe based on Category and Count
  kegg_df <- kegg_df %>% arrange(Category, Count)
  kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)
  
  # 4. Map categories to global colors
  cols_KEGG <- c(
    "Immune & Inflammation"         = uni_colors[["Red"]],
    "Signal Transduction"           = uni_colors[["Blue"]],
    "Metabolism & Cellular Process" = uni_colors[["Blue"]],
    "Pathogen Infection"            = uni_colors[["Green"]],
    "Related Pathways"              = uni_colors[["Green"]]
  )
  
  # 5. Plot generation
  pB <- ggplot(kegg_df, aes(x = Count, y = Description, fill = Category)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = Count), hjust = -0.3, size = 3) +
    scale_fill_manual(values = cols_KEGG) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) + 
    labs(title = "KEGG Enrichment Analysis", x = "Gene Count", y = NULL) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 9, color = "black"),
      
      # [Customization 2] Left-align title to the panel border
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      plot.title.position = "panel", 
      
      # Legend customization (Bottom-right inside the plot area with white background)
      legend.position = c(0.98, 0.05), 
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.9), color = "black", size = 0.2),
      legend.title = element_blank(), 
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 7),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.margin = margin(4, 4, 4, 4)
    )
  
  pdf("Figure_B_KEGG.pdf", width = 9, height = 7)
  print(pB)
  dev.off()
}

message(">>> Plotting completed! Titles are properly left-aligned to the panels.")