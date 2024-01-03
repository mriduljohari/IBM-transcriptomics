# This script is for figure 2 and 5a
# Load the workspace
load("IBMtranscriptomics.RData")

# Pacman for package management
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(RColorBrewer, gridExtra, ggplot2, ggdist, dplyr, DESeq2)

# Define genes and plot colors
genes <- c('ENSG00000006459.10;KDM7A', 'ENSG00000108679.12;LGALS3BP', 'ENSG00000164342.12;TLR3', 'ENSG00000033627.16;ATP6V0A1', 'ENSG00000140968.10;IRF8', 'ENSG00000155465.18;SLC7A7', 'ENSG00000127951.6;FGL2', 'ENSG00000166710.18;B2M', 'ENSG00000198336.9;MYL4', 'ENSG00000156587.15;UBE2L6', 'ENSG00000005206.16;SPPL2B', 'ENSG00000167552.13;TUBA1A', 'ENSG00000090447.11;TFAP4', 'ENSG00000204287.13;HLA-DRA', 'ENSG00000176887.6;SOX11')
cols_plot <- c("IBM" = "#DDAA33", "group_TMD" = "#004488", "group_Amputee" = "#CC3311")

# Common theme for presentation
blank_theme <- theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 15), legend.text = element_text(size = 15), 
        plot.title = element_text(size = 12, hjust = 0.5), axis.title.x = element_blank())

# Function to create and save plots
create_and_save_plot <- function(gene, dds, filename) {
  hgnc_symbol <- strsplit(gene, ";")[[1]][2]  # Extract HGNC symbol for display

  gene_counts <- plotCounts(dds, gene = gene, intgroup = c("Cohort"), returnData = TRUE, normalized = TRUE, transform = TRUE)
  plot <- ggplot(gene_counts, aes(x = Cohort, y = count, fill = Cohort)) +
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA, alpha = 0.35) +
    geom_boxplot(width = .2, outlier.shape = NA, alpha = 0.35) +
    geom_point(size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1)) +
    coord_cartesian(xlim = c(1.2, NA), clip = "off") +
    scale_y_log10() + scale_fill_manual(values = cols_plot) +
    blank_theme +
    ggtitle(bquote(italic(.(hgnc_symbol))))  # Corrected title expression

  ggsave(filename, plot, width = 297, height = 210, units = "mm", dpi = 300)
  return(plot)
}


# Generate plots and store in a list
plot_list <- list()
for (gene in genes) {
  plot_list[[gene]] <- create_and_save_plot(gene, dds_main, paste0("Figure_", i, ".png"))
}


# Arrange and save the panel of plots
figure <- ggarrange(plotlist = plot_list, 
                    ncol = 5, nrow = 3,
                    common.legend = TRUE, legend = "bottom",
                    font.label = list(size = 10))
annotate_figure(figure,
                left = text_grob("Normalised counts", rot = 90),
                bottom = text_grob("Cohorts", size = 10)
                )
ggsave("top15genes.png", figure, width = 297, height = 210, units = "mm", dpi = 600)

