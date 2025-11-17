######## This code compares automatic and exact cell count approahces  ##########

# GSE282765 and GSE182001
# take forced cell determination

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load objects 
obj <- readRDS("/scratch/khandl/4.Technical/Forced_cell_determination_exonic_reads_only_Hs_PB_healthy_and_CRC.rds")
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")

obj <- subset(obj, subset = percent.mt < 25)

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

#### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)                
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

##### add annotation based on cell id match 
current.cluster.ids <- c("?", "Basophils","DCs","Eosinophils","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","ProNeutro", "T")
new.cluster.ids <- c("Undefined", "Basophils","DCs","Eosinophils","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","ProNeutro", "T")
obj_reference$annotation <- plyr::mapvalues(x = obj_reference$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_reference$annotation)))$Var1

cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_reference) <- "annotation"
  sub <- subset(obj_reference, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #42240
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #1

obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Basophils ~ "Basophils",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj$annotation)
DimPlot(obj, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_reference, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

# remove NA 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Basophils","DCs","Eosinophils","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                              "Monocytes","Neutrophils","PCs","T","ProNeutro","Undefined"))

### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_exons_only_annotation.rds")



