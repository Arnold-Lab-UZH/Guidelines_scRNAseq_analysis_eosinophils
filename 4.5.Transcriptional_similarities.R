########### This code analyses transcriptional similarity of eosinophils with other cell types ##########
### Datasets used: GSM7919060, GSE216189, GSE282765, GSE182001

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Mm 
obj1 <- readRDS(file.path(seurat_objects_dir,"Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj2 <- readRDS(file.path(seurat_objects_dir,"Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj3 <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj4 <- readRDS(file.path(seurat_objects_dir,"Mm_bm_GSM7919060_anno.rds"))
obj5 <- readRDS(file.path(seurat_objects_dir,"Mm_liver_GSE216189.rds"))


obj <- merge(obj1, c(obj2,obj3,obj4,obj5))
obj <- JoinLayers(obj)

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B", "Basophils", "DCs", "EoP","Eosinophils","GMPs","HSCs", "Kupffer", "Macrophages","Mast",
                              "Monocytes","Neutrophils","PCs", "ProMono", "ProNeutro","T","TAMs"))

### Pearson correlation of HVG 
hvg <- VariableFeatures(obj)

average_expression <- AggregateExpression(obj,features = hvg, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
df <- as.data.frame(average_expression$RNA)
cor_mat <- cor(df, method = "pearson")

pheatmap(cor_mat,
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         color = colorRampPalette(c("blue", "white", "darkred"))(100),
         main = "Pearson Correlation Between Cell Types")

