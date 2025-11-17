########### This code quantifies ambient RNA content using DecontX default (cell clusters) ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

# GSE282765 
# take forced cell determination and intronic + exonic reads

##### read R objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

current.cluster.ids <- c("?","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                         "Monocytes","Neutrophils","PCs","TAMs", "T")
new.cluster.ids <- c("Undefined","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                     "Monocytes","Neutrophils","PCs","TAMs", "T")
obj$annotation <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

##### DecontX without supplying the empty droplet background, estimation of contamination based on profiels of other cell clusters in the filtered data = DecontX_default
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce <- decontX(sce,z=sce$annotation, batch = sce$experiment)

##### plot contamination in feature Plot 
# add contamination result to meta data 
obj$decontX_contamination <- colData(sce)$decontX_contamination
Idents(obj) <- "condition"
sub <- subset(obj, idents = c("P1_tumor","P2_tumor","P3_tumor"))
p <- FeaturePlot(sub, features = c("decontX_contamination")) +scale_color_gradientn( colours = c('darkblue',"yellow" ,'darkred'))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_feature_plot_Hs.svg", width = 10, height = 8, plot = p)

##### plot level of contamination based on DecontX_default 
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_decontX = median(decontX_contamination, na.rm = TRUE)) %>%
    as.data.frame()
  df$condition <- i 
  df$species <- "Mm"
  df_list[[i]] <- df
}
df <- bind_rows(df_list)

df <- df[df$celltype %in% c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                            "Monocytes","Neutrophils","PCs","TAMs", "T"),]
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/decontX_Hs.csv")

p <- ggplot(df, aes(x = reorder(celltype,median_decontX, FUN = median), y =  median_decontX, fill = celltype)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_contamination_per_celltypeHs.svg", width = 10, height = 8, plot = p)
