########## This code does sample integration, pre-processing, clustering and annotation of Hs tumor and NAT from GSE282765  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Seurat object generation 
### Forced cell determination intronic and exonic reads 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

### Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

### Add mitochondrial percentage per cell 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

### Add conditions to metadata 
patients$cell_determination <- "forced"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

### Save object
saveRDS(patients, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds")

##### Load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
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
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
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
DotPlot(obj, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" # endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:25)
new.cluster.ids <- c("Neutrophils","T","PCs","Macrophges_TAMs", "Eosinophils", "lowQ",
                     "Epithelial","Monocytes","?", "Endothelial","lowQ",
                     "Fibroblasts","PCs","Mast","Mixed","Mixed",
                     "B","T","Mast","Mixed","T",
                     "Mixed","Mixed","?","?","Mixed")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")

### Subcluster Epithelial 
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Epithelial",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.05)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Epithelial_0","Epithelial_1"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 0 = lowQ, 1 = epithelial 

# Rename
current.cluster.ids <- c("?","B","Endothelial","Eosinophils","Epithelial_0","Epithelial_1",
                         "Fibroblasts","lowQ","Macrophges_TAMs","Mast","Mixed","Monocytes", "Neutrophils","PCs","T")
new.cluster.ids <- c("?","B","Endothelial","Eosinophils","lowQ","Epithelial",
                     "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes", "Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster B cells
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "B",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.05)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "B_0","B_1"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features = "MS4A1", reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 0 = B, 1 = mixed 

# Rename
current.cluster.ids <- c("?", "B_0","B_1","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes", "Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","Mixed","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes", "Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Monocytes
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Monocytes",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Monocytes_0","Monocytes_1","Monocytes_2"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
FeaturePlot(sub_celltype, features = "CD1C")
FeaturePlot(sub_celltype, features = "VCAN")
# 2 = DCs, 0 and 1 = moncoytes

# Rename
current.cluster.ids <- c("?", "B","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes_0","Monocytes_1","Monocytes_2", "Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Monocytes","DCs", "Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Endothelial
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Endothelial",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Endothelial_0","Endothelial_1","Endothelial_2"))
DimPlot(sub_celltype, reduction = "umap.mnn")

FeaturePlot(sub_celltype, features =c("PECAM1","COL1A2"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 1 = Endothelial, 0 and 2 = Fibroblasts

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial_0","Endothelial_1","Endothelial_2","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","DCs", "Fibroblasts","Fibroblasts","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Fibroblasts
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Fibroblasts",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Fibroblasts_0","Fibroblasts_1","Fibroblasts_2","Fibroblasts_3"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features =c("PECAM1","COL1A2"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
#3 = Mixed, 2 = Endothlelial; 1 and 0 = Fibroblasts

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts_0","Fibroblasts_1","Fibroblasts_2","Fibroblasts_3",
                        "lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts","Fibroblasts","Endothelial","Mixed",
                     "lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Eosinophils
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Eosinophils",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Eosinophils_0","Eosinophils_1","Eosinophils_2"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features = c("CLC","KIT","CLEC12A","TPSB2","SYNE1","FFAR2"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 0 = Eosinophils, rest is lowQ

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils_0","Eosinophils_1","Eosinophils_2","Epithelial","Fibroblasts",
                         "lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","lowQ","lowQ","Epithelial","Fibroblasts",
                     "lowQ","Macrophages_TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Macrophages_TAMs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Macrophages_TAMs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Macrophages_TAMs_0","Macrophages_TAMs_1"))
DimPlot(sub_celltype, reduction = "umap.mnn")

FeaturePlot(sub_celltype, features = c("C1QC","FOLR2","SPP1"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 1 = TAMs, 2 = Macrophages 

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages_TAMs_0","Macrophages_TAMs_1","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","TAMs","Mast","Mixed","Monocytes","Neutrophils","PCs","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Mast
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mast",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.05)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mast_0","Mast_1","Mast_2", "Mast_3"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features = c("CPA3","KIT","CLEC12A","CD200R1","KLK10"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# Mast only 0 , everything else mixed 

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast_0","Mast_1","Mast_2", "Mast_3","Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Mixed", "Mixed","Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "PCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "PCs_0","PCs_1","PCs_2","PCs_3"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features = c("IGKC","IGHA2","JCHAIN"))
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 0,1 = PCs, everything else is mixed 

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs_0","PCs_1","PCs_2","PCs_3","T","TAMs")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","PCs","Mixed","Mixed","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster T cells
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_0","T_1","T_2","T_3"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
FeaturePlot(sub_celltype, features = c("CD3E","CD4","CD8A","ICOS","IL17A","GZMA"), reduction = "umap.mnn")
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
# 0 = T, rest is mixed 

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T_0","T_1","T_2","T_3","TAMs")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","Mixed","Mixed","Mixed","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### Check annotation
Idents(subCl) <- "annotation"

DotPlot(subCl, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1","FFAR2", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts
                                 "PECAM1" # Endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
obj <- subCl
obj$sample <- obj$condition
obj$tissue2 <- obj$tissue
current.cluster.ids <- c("tissue_ctrl","tumor")
new.cluster.ids <- c("colon","tumor")
obj$tissue <- plyr::mapvalues(x = obj$tissue2, from = current.cluster.ids, to = new.cluster.ids)
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

