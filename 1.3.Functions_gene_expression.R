##### Functions to analyse gene expreession (DEGs, correlation)
### Function to do DGE analysis between two conditions and writing the output to a csv file 
DEG_to_csv_two_cond <- function(
    seurat_object,
    assay_oi,
    cond1,
    cond2,
    only.pos_condition,
    logfc_threshold,
    csv_file_directory
){
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = assay_oi)
  DefaultAssay(seurat_object) <- assay_oi
  markers <- FindMarkers(object = seurat_object, ident.1 = cond1, ident.2 = cond2, only.pos = only.pos_condition, min.pct = 0.25, 
                         logfc.threshold = logfc_threshold,slot = "data")
  write.csv(markers, file = csv_file_directory)
}

DEG_two_cond_pb_DESeq2 <- function(
    pseudobulk_object,
    celltype_cond1,
    celltype_cond2,
    output_path
){
  bulk_de <- FindMarkers(object = pseudobulk_object, 
                         ident.1 = celltype_cond1,
                         ident.2 = celltype_cond2,
                         test.use = "DESeq2")
  # write output 
  write.csv(bulk_de, output_path)
}

### Function to plot genes of interest per condition of interest, scaled per row 
heatmap_goi_coi <- function(
    seurat_object,
    condition_oi,
    markers_oi,
    groups_of_markers,
    number_of_markers_per_group,
    colors_per_group,
    groups_and_colors,
    cluster_rows_cond,
    cluster_cols_cond
){
  # average expression per cluster and condition 
  average_expression <- AverageExpression(seurat_object, return.seurat = FALSE, features = markers_oi, normalization.method = "LogNormalize",assays = "RNA", group.by = condition_oi)
  average_expression_df <- as.data.frame(average_expression)
  average_expression_df <- average_expression_df[match(markers_oi, rownames(average_expression_df)),]
  
  # prepare palette for pheatmap
  paletteLength   <- 50
  myColor         <- colorRampPalette(c("darkblue", "white","orange" ,"darkred"))(paletteLength)
  breaksList      = seq(-2, 2, by = 0.04)
  
  # prepare annotation for pheatmap
  annotation_rows             <- data.frame(markers = rep(groups_of_markers, number_of_markers_per_group))
  rownames(annotation_rows)   <- rownames(average_expression_df)
  annotation_rows$markers     <- factor(annotation_rows$markers, levels = groups_of_markers)
  
  mycolors <- colors_per_group
  names(mycolors) <- unique(annotation_rows$markers)
  mycolors <- list(category = mycolors)
  annot_colors=list(markers=groups_and_colors)
  
  p <- pheatmap(average_expression_df,scale = "row",
                color = colorRampPalette(c("darkblue", "white","orange", "darkred"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                breaks = breaksList,
                cluster_rows = cluster_rows_cond, cluster_cols = cluster_cols_cond, 
                border_color = "black", 
                legend_breaks = -6:6, 
                cellwidth = 10, cellheight = 5,
                angle_col = "45", 
                annotation_colors = annot_colors,
                annotation_row = annotation_rows,
                fontsize = 5)
  print(p)
}
