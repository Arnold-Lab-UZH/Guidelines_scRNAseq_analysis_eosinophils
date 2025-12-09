########## This code uses DecontX for ambient RNA quantification  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load R object 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Run DecontX 
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce <- decontX(sce,z=sce$annotation, batch = sce$experiment)

##### plot contamination in feature Plot 
### Add the contamination result to the meta.data from the Seurat object 
obj$decontX_contamination <- colData(sce)$decontX_contamination
Idents(obj) <- "condition"
sub <- subset(obj, idents = c("P1_tumor","P2_tumor","P3_tumor"))
p <- FeaturePlot(sub, features = c("decontX_contamination")) +scale_color_gradientn( colours = c('darkblue',"yellow" ,'darkred'))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_feature_plot_Hs.svg", width = 10, height = 8, plot = p)

### Generate a dataframe of decontX_contamination for each sample and cell types 
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
### save data frame 
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/decontX_Hs.csv")

