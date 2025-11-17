########## This code quantifies ambient RNA content  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

# GSE282765 and GSE182001
# take forced cell determination and intronic + exonic reads

##### read R objects 
human_colon_tumor <- readRDS( "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
mouse_colon_tumor <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

##### plot IGKC and IGLC2 in feature plot 
Idents(human_colon_tumor) <- "condition"
sub <- subset(human_colon_tumor, idents = c("P1_tumor","P2_tumor","P3_tumor"))
p <- FeaturePlot(sub, features = c("IGKC")) +scale_color_gradientn( colours = c('grey', 'darkred'))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/IGKC_human.svg", width = 10, height = 8, plot = p)

Idents(mouse_colon_tumor) <- "condition"
sub <- subset(mouse_colon_tumor, idents = "tumor_wt")
p <- FeaturePlot(sub, features = c("Igkc")) +scale_color_gradientn( colours = c('grey', 'darkred'),limits = c(0,8.5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/IGKC_mouse.svg", width = 10, height = 8, plot = p)

##### Average expression of Igkc in Boxplot 
### human_colon 
obj <- human_colon_tumor
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  average_expression <- AverageExpression(sub, features = "IGKC", return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
  average_expression_df <- as.data.frame(average_expression)
  average_expression_df <- as.data.frame(t(average_expression_df))
  average_expression_df$condition <- i 
  average_expression_df$species <- "Hs"
  average_expression_df$celltype <- rownames(average_expression_df)
  colnames(average_expression_df) <- c("average_expr","condition","species", "celltype")
  df_list[[i]] <- average_expression_df
}
df1 <- bind_rows(df_list)

### mouse_colon
obj <- mouse_colon_tumor
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  average_expression <- AverageExpression(sub, features = "Igkc", return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
  average_expression_df <- as.data.frame(average_expression)
  average_expression_df <- as.data.frame(t(average_expression_df))
  average_expression_df$condition <- i 
  average_expression_df$species <- "Mm"
  average_expression_df$celltype <- rownames(average_expression_df)
  colnames(average_expression_df) <- c("average_expr","condition","species", "celltype")
  df_list[[i]] <- average_expression_df
}
df2 <- bind_rows(df_list)

df_list_all <- list(df1,df2)
df <- bind_rows(df_list_all)

df <- df[df$celltype %in% c("RNA.Eosinophils","RNA.Epithelial","RNA.Neutrophils","RNA.Fibroblasts","RNA.B","RNA.T","RNA.Monocytes",
                            "RNA.Mast","RNA.PCs","RNA.Endothelial","RNA.Macrophages","RNA.TAMs","RNA.DCs"),]
p <- ggplot(df, aes(x = reorder(celltype,average_expr, FUN = median), y =  average_expr, fill = celltype)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + scale_y_break(c(200, 2000))   + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values = c("RNA.Eosinophils" = "#E22F27", "RNA.Epithelial" = "#6D5421",  "RNA.Neutrophils" = "#9518ED",  "RNA.Fibroblasts" = "#443511", 
                               "RNA.B" = "#F3E972", "RNA.T" = "#5BC7D9",   "RNA.Monocytes" = "#ADD8AB", "RNA.Mast" = "#7F7F79", "RNA.PCs" = "#B4C108",
                               "RNA.Endothelial" = "#A09167", "RNA.Macrophages" = "#82C341", 
                               "RNA.TAMs" = "#516D38","RNA.DCs" = "#E43794","RNA.lowQ"="#282525","RNA.Mixed" = "#E59A38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/IGKC_averag_expression_all.svg", width = 10, height = 8, plot = p)

## statistical test --> one way ANOVA 
anova <- aov(average_expr ~ celltype, data = df)
summary(anova)
TukeyHSD(anova)
