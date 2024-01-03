# This script is for figure 1
# load the workspace if not already loaded
# load.image('IBMtranscriptomics.RData')
# Ensure all required packages are installed and loaded
if (!require(pacman)) install.packages("pacman")
pacman::p_load(RColorBrewer, gridExtra, ggplot2, VennDiagram, dplyr, ggbeeswarm)

# Define a function to save plots
save_plot <- function(plot, filename, width, height, pointsize, res) {
  png(filename, width = width, height = height, pointsize = pointsize, res = res, bg = "white")
  grid.arrange(plot, ncol = 2, top = textGrob("Principal component analysis", gp = gpar(fontsize = 15, font = 8)))
  dev.off()
}

# PCA plot - Assuming p1 and p2 are already created PCA plots
save_plot(p1, "Principal_component_analysis_IBM_AMP.png", 4800, 3600, 8, 600)
save_plot(p2, "Principal_component_analysis_IBM_TMD.png", 4800, 3600, 8, 600)  # If p2 is a separate plot

# Define a function for Venn diagrams for both read.csv and readLines inputs
create_venn <- function(genes1, genes2, label1, label2, filename, colors) {
  venn.diagram(
    x = list(label1 = genes1, label2 = genes2),
    category.names = c(label1, label2),
    filename = filename,
    output = TRUE,
    imagetype = "png",
    units = 'px',
    height = 4800,
    width = 4800,
    resolution = 600,
    lty = 'blank',
    fill = colors,
    margin = 0.07,
    col = "transparent",
    alpha = c(0.3, 0.3),
    cex = 2.2,
    euler.d = TRUE,
    scaled = TRUE,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.pos = c(335, -195),
    cat.fontface = "bold",
    cat.fontfamily = "sans"
  )
}


# Venn diagrams for up and downregulated genes
# Read files and extract first column
up_IBM_AMP_genes <- read.csv("resSort_main_up.lfcShrink.IBM_AMP.csv")[,1]
up_TMD_AMP_genes <- read.csv("resSort_main_up.lfcShrink.TMD_AMP.csv")[,1]
down_IBM_AMP_genes <- read.csv("resSort_main_down.lfcShrink.IBM_AMP.csv")[,1]
down_TMD_AMP_genes <- read.csv("resSort_main_down.lfcShrink.TMD_AMP.csv")[,1]

# Create Venn diagrams
create_venn(up_IBM_AMP_genes, up_TMD_AMP_genes, "IBMvAmputee", "TMDvAmputee", "Main_IBMspecific_upregulated_DEseq2_final.png", c('yellow', 'green'))
create_venn(down_IBM_AMP_genes, down_TMD_AMP_genes, "IBMvAmputee", "TMDvAmputee", "Main_IBMspecific_downregulated_DEseq2.png", c('yellow', 'green'))

# Venn diagram for DGE and DSG lists
DGE_list <- readLines("DGE_list_for Venn.txt")
DSG_list <- readLines("DSG_list_for Venn.txt")
create_venn(DGE_list, DSG_list, "IBM-specific differentially expressed genes", "IBM-specific differentially spliced genes", "DGE_DSG_genes_venn.png", c('cyan', 'magenta'))


# Save workspace for later reuse
save.image('IBMtranscriptomics.RData')