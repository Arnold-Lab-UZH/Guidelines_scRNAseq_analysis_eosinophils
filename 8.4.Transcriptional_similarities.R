########## This code analyzes transcriptional correlation between cell types  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Mm 
obj1 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
obj2 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
obj3 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")
obj4 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_bm_GSM7919060_anno.rds")
obj5 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_liver_GSE216189.rds")

obj <- merge(obj1, c(obj2,obj3,obj4,obj5))
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B", "Basophils", "DCs", "EoP","Eosinophils","GMPs","HSCs", "Kupffer", "Macrophages","Mast",
                              "Monocytes","Neutrophils","PCs", "ProMono", "ProNeutro","T","TAMs"))


### Pearson correlation of HVG 
hvg <- FindVariableFeatures(obj)
hvg <- VariableFeatures(obj)

average_expression <- AggregateExpression(obj,features = hvg, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
df <- as.data.frame(average_expression)
df <- data.frame(lapply(df, as.numeric), row.names = rownames(df))

cor_mat <- cor(df, method = "pearson")

pheatmap(cor_mat,
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         color = colorRampPalette(c("blue", "white", "darkred"))(100),
         main = "Pearson Correlation Between Cell Types")

