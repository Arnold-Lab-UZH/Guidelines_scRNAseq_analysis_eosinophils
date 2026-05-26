########## This code does sample integration, pre-processing, clustering and annotation of Mm CRC tumor, NAT and disseminated and healthy colon from GSE282765  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Seurat object generation 
### forced cell determination - intronic and exonic reads 
tumor_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_tumor_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_disseminated_ST07_Expression_Data.st"), 
  project = "disseminated_wt", condition = "disseminated_wt",3,200)

adult_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_colon_ST04_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

adjacent_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_NAT_ST03_Expression_Data.st"), 
  project = "adjacent_colon_wt", condition = "adjacent_colon_wt",3,200)

### Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt,adjacent_colon_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt","adjacent_colon_wt"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

### Add conditions to metadata 
tumor$cell_determination <- "forced"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

### Save object
saveRDS(tumor, file = "/scratch/khandl/4.Technical/Forced_cell_determination_intronic_and_exonic_reads_Mm_CRC_AKPS_tumor_NAT_diss.rds")

##### Load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_intronic_and_exonic_reads_Mm_CRC_AKPS_tumor_NAT_diss.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)                     
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### nFeature and percent.mito per cluster to exclude low quality clusters 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

### Marker gene expression 
DotPlot(obj, features =  c("Tpsab1", "Tpsb2","Kit",# Mast cells 
                           "Mcpt8","Cd200r3","Clec12a", #B asophils
                           "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                           "Clec9a","Flt3","Xcr1" ,# DCs
                           "Cd300e","Ereg","Vcan", # Monocytes
                           "F13a1", "Folr2","C1qc" ,# Macrophages 
                           "Spp1","Mmp12","Fn1", # TAMs
                           "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                           "Icos","Cd8a", "Cd3e","Trac", # T cells 
                           "Cd19", "Cr2","Ms4a1", # B cells 
                           "Igha","Jchain","Igkc", # PCs
                           "S100a8","S100a9", #Neutrophils 
                           "Epcam","Tff3","Muc2", # Epithelial 
                           "Col1a2", # Fibroblasts/
                           "Pecam1" # Endothelial 
)) + theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:21)
new.cluster.ids <- c("Neutrophils","T","TAMs","lowQ", "B", "Neutrophils",
                     "Macrophages","Eosinophils","PCs","Baso_Mast_mixed","lowQ",
                     "?","DCs","Mixed","PCs","Mixed",
                     "Mixed","Mixed","Mixed","?","?","Mixed")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Baso_Mast_mixed
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Baso_Mast_mixed",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.8)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Baso_Mast_mixed_0","Baso_Mast_mixed_1","Baso_Mast_mixed_2",
                                         "Baso_Mast_mixed_3","Baso_Mast_mixed_4","Baso_Mast_mixed_5", "Baso_Mast_mixed_6"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# Rename
current.cluster.ids <- c( "?","B","Baso_Mast_mixed_0","Baso_Mast_mixed_1","Baso_Mast_mixed_2",
                          "Baso_Mast_mixed_3","Baso_Mast_mixed_4","Baso_Mast_mixed_5", "Baso_Mast_mixed_6",
                          "DCs", "Eosinophils","lowQ",
                          "Macrophages","Mixed","Neutrophils","PCs","T","TAMs")
new.cluster.ids <- c(  "?","B","Mixed","Mixed","Mixed",
                       "Mixed","Mixed","Mixed", "Mast",
                       "DCs", "Eosinophils","lowQ",
                       "Macrophages","Mixed","Neutrophils","PCs","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster TAMs to find monoytes 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "TAMs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "TAMs_0","TAMs_1"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 1 = monocytes, 0 = TAMs

# Rename
current.cluster.ids <- c(  "?","B", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Neutrophils","PCs","T", "TAMs_0","TAMs_1")
new.cluster.ids <- c( "?","B", "DCs", "Eosinophils","lowQ",
                      "Macrophages","Mast", "Mixed","Neutrophils","PCs","T", "TAMs","Monocytes")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### subcluster T cells
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_0","T_1","T_2","T_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 3 = Mixed, all other T 

# Rename
current.cluster.ids <- c(  "?","B", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes", "Neutrophils","PCs","T_0","T_1","T_2","T_3","TAMs")
new.cluster.ids <- c(   "?","B", "DCs", "Eosinophils","lowQ",
                        "Macrophages","Mast", "Mixed","Monocytes", "Neutrophils","PCs","T","T","T","Mixed","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster B cells
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "B",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "B_0","B_1","B_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 2 = lowQ, all other B 

# Rename
current.cluster.ids <- c(  "?", "B_0","B_1","B_2", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
new.cluster.ids <- c(   "?", "B","B","lowQ", "DCs", "Eosinophils","lowQ",
                        "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "PCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "PCs_0","PCs_1","PCs_2","PCs_3","PCs_4","PCs_5"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
VlnPlot(sub_celltype, features = "nFeature_RNA")
FeaturePlot(sub_celltype, features = "Igkc")
# 0,1,2,4 = PCs, the rest is Mixed 

# Rename
current.cluster.ids <- c(  "?", "B","DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs_0","PCs_1","PCs_2","PCs_3","PCs_4","PCs_5","T","TAMs")
new.cluster.ids <- c(  "?", "B","DCs", "Eosinophils","lowQ",
                       "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","PCs","PCs","Mixed","PCs","Mixed","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Check annotation
obj <- subCl
Idents(obj) <- "annotation"
DotPlot(obj, features =  c("Tpsab1", "Tpsb2","Kit", # Mast cells
                           "Mcpt8","Cd200r3","Clec12a", # Basophils
                           "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
                           "Clec9a","Flt3","Xcr1" ,# DCs
                           "Cd300e","Ereg","Vcan", # Monocytes
                           "F13a1", "Folr2","C1qc" ,# Macrophages 
                           "Spp1","Mmp12","Fn1", # TAMs
                           "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                           "Icos","Cd8a", "Cd3e","Trac", # T cells 
                           "Cd19", "Cr2","Ms4a1", # B cells 
                           "Jchain","Igkc", # PCs
                           "S100a8","S100a9", # Neutrophils 
                           "Epcam","Tff3","Muc2", # Epithelial 
                           "Col1a2", # Fibroblasts/
                           "Pecam1" # Endothelial 
)) + theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
obj$sample <- obj$condition
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids <- c("colon","colon","tumor","tumor")
obj$tissue <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
