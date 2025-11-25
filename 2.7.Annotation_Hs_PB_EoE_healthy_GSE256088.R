########## This code does sample integration, pre-processing, clustering and annotation of Hs PB from healthy and EoE patients from GSE256088  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Seurat object generation 
PB_ctrl1 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_ctrl1","EoE",3,200, "PB_ctrl1","PB_healhty","EoE_Exp1","healthy")
PB_ctrl2 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_ctrl2","EoE",3,200, "PB_ctrl2","PB_healhty","EoE_Exp2","healthy")
PB_ctrl3 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_ctrl3","EoE",3,200, "PB_ctrl3","PB_healhty","EoE_Exp3","healthy")

PB_EoE1 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_EoE1","EoE",3,200, "PB_EoE1","PB_EoE","EoE_Exp1","EoE")
PB_EoE2 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_EoE2","EoE",3,200, "PB_EoE2","PB_EoE","EoE_Exp2","EoE")
PB_EoE3 <- create_seurat_10X_structured_data( "/scratch/khandl/2.EoE_Morgenstern2024_Morgan2021/GSE256088_RAW/PB_EoE3","EoE",3,200, "PB_EoE3","PB_EoE","EoE_Exp3","EoE")

### Merge samples
blood <- merge(PB_ctrl1, y = c(PB_ctrl2, PB_ctrl3,PB_EoE1,PB_EoE2,PB_EoE3),
               add.cell.ids = c("PB_ctrl1", "PB_ctrl2","PB_ctrl3","PB_EoE1", "PB_EoE2","PB_EoE3"))
blood <- JoinLayers(blood)

### Add mitochondrial percentage per cell 
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

### Apply mitochondrial cutoff 
blood <- subset(blood, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- blood
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)                     
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap.mnn",group.by = "condition")
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

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
                                 "CLC","ADGRE1","SYNE1","FFAR2", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" , #Endothelial 
                                 "CD34","CEBPA","MSI1", # GMPs
                                 "EPX","EAR2","EAR1", # EoP
                                 "ELANE","CEBPE","S100A8" ,# ProNeutro
                                 "LY6C2","CSF1R","IRF8" # ProMono
)) )+ theme(axis.text.x = element_text(angle = 90)) 

FeaturePlot(obj, features = c("CLC","VCAN","CD3E","FCGR3B","CD79A"), reduction = "umap.mnn")

### Rename clusters
current.cluster.ids <- c(0:6)
new.cluster.ids <- c("lowQ","lowQ","Neutrophils","T", "Monocytes","Eosinophils","B")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,reduction = "umap.mnn" )

##### Check annotation
Idents(obj) <- "annotation"
DotPlot(obj, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
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
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" , #Endothelial 
                                 "CD34","CEBPA","MSI1", # GMPs
                                 "EPX","EAR2","EAR1", # EoP
                                 "ELANE","CEBPE","S100A8" ,# ProNeutro
                                 "LY6C2","CSF1R","IRF8" # ProMono
)) )+ theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
obj$technology <- "HIVE"
obj$species <- "Hs"
obj$cell_enrichment <- "Eosinophils"
obj$tissue <- "blood"

##### Save object 
saveRDS(obj, file = "/scratch/khandl/technical/seurat_objects/Hs_blood_GSE256088_et_al_anno.rds")

