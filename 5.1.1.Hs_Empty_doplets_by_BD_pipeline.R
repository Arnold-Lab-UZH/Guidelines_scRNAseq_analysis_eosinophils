########### This code compares forced and automatic BD pipeline ##########

# Hs tumor and NAT GSE282765
# take forced and automatic cell determination and intronic + exonic reads

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load objects 
obj_forced <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
obj_autom <- readRDS( "/scratch/khandl/4.Technical/Automatic_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds")

##### pre-process and cluster 
obj <- obj_autom
obj <- subset(obj, subset = percent.mt < 25)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)         
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)
obj_autom <- obj

##### add annotation based on cell id match 
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
obj_forced$annotation <- plyr::mapvalues(x = obj_forced$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_forced$annotation)))$Var1

cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_forced) <- "annotation"
  sub <- subset(obj_forced, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

shared <- intersect(rownames(obj_autom@meta.data),  rownames(obj_forced@meta.data)) #78669
additional_ones_in_autom <- rownames(obj_autom@meta.data)[!rownames(obj_autom@meta.data) %in% shared] #673

obj_autom$shared_unique <- NA
obj_autom@meta.data <- obj_autom@meta.data %>%
  mutate(shared_unique = case_when(
    rownames(obj_autom@meta.data) %in% shared ~ "shared_BC",
    rownames(obj_autom@meta.data)  %in% additional_ones_in_autom ~ "unique_BD_autom",
    TRUE ~ NA_character_))
table(obj_autom$shared_unique)
DimPlot(obj_autom, group.by = "shared_unique")

obj_autom$annotation <- NA
obj_autom@meta.data <- obj_autom@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(obj_autom$annotation)
DimPlot(obj_autom, group.by = "annotation", label = TRUE)

## susbset additional_ones_in_autom and annotate
Idents(obj_autom) <- "shared_unique"
sub <- subset(obj_autom, idents = "unique_BD_autom")
DimPlot(sub, group.by = "mnn.clusters", label =TRUE)

Idents(sub) <- "mnn.clusters"
DotPlot(sub, features = unique(c("TPSAB1","KIT","TPSB2", #Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,#DCs
                                 "CD300E","EREG","VCAN", #Monocytes
                                 "F13A1", "FOLR2","C1QC" ,#Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1","FFAR2", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", #Neutrophils 
                                 "EPCAM","TFF3","MUC2", #Epithelial 
                                 "COL1A2","PECAM1" # Fibroblasts/endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 
VlnPlot(sub, features = "nFeature_RNA")

# rename
current.cluster.ids <- c(0:28)
new.cluster.ids <- c("PCs","T","Neutrophils", "Macrophages","lowQ", "TAMs", 
                     "Monocytes","Epithelial","Fibroblasts","PCs", "Fibroblasts",
                     "Mixed","Mast","T","Mixed","Mixed",
                     "Mixed","Neutrophils","Mixed","Mixed","Mixed",
                     "Mixed","Mixed","Mixed","Mast","Mixed",
                     "Mixed","Mixed","Mixed")
sub$annotation <- plyr::mapvalues(x = sub$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
p1 <- DimPlot(sub, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj_autom, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

## merge objects 
# remove NA from obj_autom
Idents(obj_autom) <- "annotation"
obj_autom <- subset(obj_autom, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                                          "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

obj_autom <- merge(obj_autom, sub)
obj_autom <- JoinLayers(obj_autom)

### save object 
saveRDS(obj_autom, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_automatic_cell_determination_with_intronic_reads_annotated.rds")



