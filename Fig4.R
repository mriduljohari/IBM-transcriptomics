# This script is for figure 4
# Load the workspace
load("IBMtranscriptomics.RData")

# Ensure all required packages are installed and loaded using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(RColorBrewer, gridExtra, ggplot2, clusterProfiler, org.Hs.eg.db, reshape2, dplyr, ComplexUpset)

# Upset plot preparation
filelist <- list.files(pattern = "*list.txt")
res <- lapply(filelist, function(x) {
  data.frame(set = x, geneID = as.character(read.table(x)[,1]), val = 1)
})
res <- dplyr::bind_rows(res)

# Rename columns
rename_map <- c(
  'Calcium_TLA_list.txt'  = 'Calcium-induced T Lymphocyte Apoptosis',
  'DSG_list.txt'          = 'IBM specific differentially spliced',
  'flux_list.txt'         = 'Flux of Ca2+',
  'mobilisation_list.txt' = 'Mobilisation of Ca2+',
  'quantity_list.txt'     = 'Quantity of Ca2+',
  'release_list.txt'      = 'Release of Ca2+'
)
# Bind rows and rename 'set' column
res <- dplyr::bind_rows(res)
res$set <- rename_map[res$set]

# Cast to wide format
res_wide <- reshape2::dcast(res, geneID ~ set, value.var = "val", fill = 0)
res_wide$name <- rownames(res_wide)

# Prepare columns for UpSet plot
upset_cols <- unname(rename_map)

# Upset plot with additional arguments
upset_plot <- ComplexUpset::upset(
res_wide, upset_cols, name = "Sets", 
    width_ratio = 0.18, height_ratio = 0.5, 
    sort_sets = 'ascending',
    set_sizes=upset_set_size(
        geom = geom_bar(width=0.2, fill ='#882255'),
        filter_intersections=FALSE,
        position='left') +
        geom_text(aes(label=..count..), hjust=1.1, stat='count') +
        expand_limits(y=1750),
    wrap=TRUE,
    base_annotations=list(
        'Intersection size'=intersection_size(
            text=list(size=2.8))),
    themes=upset_modify_themes(list(
        'Intersection size'=theme(axis.text.y=element_text(size=14),
                                  axis.title.y=element_text(size=14)),
        'intersections_matrix'=theme(text=element_text(size=16)),
        'overall_sizes'=theme(axis.text.x=element_text(size=14,angle=90))
    ))
)

ggsave("upset_plot.png", upset_plot, width = 297, height = 210, units = "mm", dpi = 300)


# Function to process gene lists for enrichGO
process_gene_list_for_enrichGO <- function(file_path) {
  gene_data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  gene_ids <- as.character(gene_data$Ensembl_id)
  return(gene_ids)
}

# Load gene lists for enrichGO
geneList_DSE <- process_gene_list_for_enrichGO("input_clusterP_IBM_DSE_ordered.txt")
genebackground_list <- as.character(read.table("input_clusterP_background.txt", header = TRUE, stringsAsFactors = FALSE)$Ensembl_id)


# GO enrichment analysis
enrich_analysis <- function(gene_list, gene_background, ont) {
  enrichGO(gene = gene_list, 
           keyType = "ENSEMBL", 
           universe = gene_background, 
           OrgDb = org.Hs.eg.db, 
           ont = ont, 
           minGSSize = 10, 
           maxGSSize = 1000, 
           pAdjustMethod = "BH", 
           pvalueCutoff = 0.01, 
           qvalueCutoff = 0.05, 
           readable = TRUE)
}

# Run enrichment analysis
ego_BP <- enrich_analysis(geneList_DSE, genebackground_list, "BP")
ego_CC <- enrich_analysis(geneList_DSE, genebackground_list, "CC")
ego_MF <- enrich_analysis(geneList_DSE, genebackground_list, "MF")

# Dot plot function
dot_plot <- function(ego_data, title) {
  dotplot(ego_data, showCategory = 20, x = "GeneRatio") +
    labs(title = "Enrichment of Gene Ontology terms", subtitle = title) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
}

# Generate and save dot plots
ggsave("dotplot_cluster_BP.png", dot_plot(ego_BP, "Biological Processes"), width = 297, height = 210, units = "mm", dpi = 300)
ggsave("dotplot_cluster_CC.png", dot_plot(ego_CC, "Cellular Component"), width = 297, height = 210, units = "mm", dpi = 300)
ggsave("dotplot_cluster_MF.png", dot_plot(ego_MF, "Molecular Function"), width = 297, height = 210, units = "mm", dpi = 300)

# save workspace for later reuse
save.image("IBMtranscriptomics.RData")