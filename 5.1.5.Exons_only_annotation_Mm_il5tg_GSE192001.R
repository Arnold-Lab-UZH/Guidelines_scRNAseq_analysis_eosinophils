########### This code compares forced BD pipeline with two different gene mapping strategies (exons + introns and exons only )  ##########
### Datasets used: GSE182001

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Load annotated object from BD forced pipeline with intron and exon mapping 
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Seurat object generation from BD forced exons only  
### Forced cell determination exonic reads only 
### forced cell determination - exonic reads only 
stomach <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_stomach_ST01_Expression_Data.st"), 
  project = "steady_state", condition = "stomach",3,200)

colon <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_colon_ST02_Expression_Data.st"), 
  project = "steady_state", condition = "colon",3,200)

small_int <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_SI_ST03_Expression_Data.st"), 
  project = "steady_state", condition = "small_int",3,200)

spleen <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_spleen_ST04_Expression_Data.st"), 
  project = "steady_state", condition = "spleen",3,200)

blood <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_blood_ST08_Expression_Data.st"), 
  project = "steady_state", condition = "blood",3,200)

bm <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_bm_ST06_Expression_Data.st"), 
  project = "steady_state", condition = "bm",3,200)

### Merge samples
merged <- merge(stomach, y = c(colon,small_int,spleen,  blood, bm),
                add.cell.ids = c("stomach","colon","SI","spleen", "blood","bm"))
meged <- JoinLayers(merged)

### Add mitochondrial percentage per cell  
meged$percent.mt <- PercentageFeatureSet(meged, pattern = "^mt.")

### Add conditions to metadata 
meged$cell_determination <- "forced"
meged$reads <- "exonic_only"
meged$species <- "Mm"
meged$technology <- "BD_Rhapsody"
meged$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(meged, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Mm_il5tg_steady_state.rds")

##### Load exons only object 
obj <- readRDS("/scratch/khandl/4.Technical/Forced_cell_determination_exonic_reads_only_Mm_il5tg_steady_state.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15)
DimPlot(obj,reduction = "umap",group.by = "seurat_clusters",raster=FALSE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

##### Transfer of annotation based on matching cell IDs 
### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?", "Endothelial","EoP", "Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mixed","Neutrophils", "PCs","ProNeutro")
new.cluster.ids <- c("Undefined","Endothelial","EoP", "Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mixed","Neutrophils", "PCs","ProNeutro")
obj_reference$annotation <- plyr::mapvalues(x = obj_reference$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_reference$annotation)))$Var1

## Loop through cell type labels and extract cell IDs for each label 
cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_reference) <- "annotation"
  sub <- subset(obj_reference, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

## Check if the cell type labels are matching between analyses 
shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #12328
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #0

## Add cell type labels based on cell IDs from BD forced 
obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$EoP ~ "EoP",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj$annotation)
DimPlot(obj, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_reference, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap")
p2 <- DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap")
p1 + p2

# remove NA 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Undefined","Endothelial","EoP", "Eosinophils","Epithelial","Fibroblasts",
                              "lowQ", "Macrophages","Mixed","Neutrophils", "PCs","ProNeutro"))

### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_exons_only_annotation.rds")

########## Annotate eosinophil subtypes additionally ##########
### Extract only Eosinophils and EoPs from newly annotated object 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Eosinophils","EoP"))

### Load the reference (intron and exon mapping )
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_annotation.rds")

##### Clustering 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### FastMNN integraiton 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:20, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 10)

##### Transfer of annotation based on matching cell IDs 
cell_types <- c("circulating","active","basal","immature","progenitor")

## Loop through cell type labels and extract cell IDs for each label 
cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_reference) <- "annotation"
  sub <- subset(obj_reference, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

## Check if the cell type labels are matching between analyses 
shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #7310
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #1264

## Add cell type labels based on cell IDs from BD forced 
obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$circulating ~ "circulating",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$active ~ "active",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$basal ~ "basal",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$immature ~ "immature",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$progenitor ~ "progenitor",
    TRUE ~ NA_character_))
table(obj$annotation)
DimPlot(obj, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_reference, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

# Remove NA
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("circulating","active","basal","immature","progenitor"))

##### Save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_forced_cell_determination_exons_only_annotation.rds")


