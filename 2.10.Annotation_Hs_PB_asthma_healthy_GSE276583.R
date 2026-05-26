########## This code does sample integration, pre-processing, clustering and annotation of PB from healthy individuals and asthma patients from GSE276583  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Seurat object generation 
PB_ctrl1 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_healthy1","healthy",3,200, "healthy_PB1","healhty_PB","Exp1","healthy")
PB_ctrl2 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_healthy2","healthy",3,200, "healthy_PB2","healhty_PB","Exp1","healthy")
PB_ctrl3 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_healthy3","healthy",3,200, "healthy_PB3","healhty_PB","Exp1","healthy")
PB_ctrl4 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_healthy4","healthy",3,200, "healthy_PB4","healhty_PB","Exp1","healthy")

PB_asthma_mild1 <-  create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma1","mild_asthma",3,200, "mild_asthma_PB1","mild_asthma","Exp1","mild_asthma")
PB_asthma_mild2 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma2","mild_asthma",3,200, "mild_asthma_PB2","mild_asthma","Exp1","mild_asthma")
PB_asthma_mild3 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma3","mild_asthma",3,200, "mild_asthma_PB3","mild_asthma","Exp1","mild_asthma")
PB_asthma_mild4 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma4","mild_asthma",3,200, "mild_asthma_PB4","mild_asthma","Exp1","mild_asthma")
PB_asthma_mild5 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma5","mild_asthma",3,200, "mild_asthma_PB5","mild_asthma","Exp1","mild_asthma")
PB_asthma_mild6 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_mild_asthma6","mild_asthma",3,200, "mild_asthma_PB6","mild_asthma","Exp1","mild_asthma")

PB_asthma_severe1 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_severe_asthma1","severe_asthma",3,200, "severe_asthma_PB1","severe_asthma","Exp1","severe_asthma")
PB_asthma_severe2 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_severe_asthma2","severe_asthma",3,200, "severe_asthma_PB2","severe_asthma","Exp1","severe_asthma")
PB_asthma_severe3 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_severe_asthma3","severe_asthma",3,200, "severe_asthma_PB3","severe_asthma","Exp1","severe_asthma")
PB_asthma_severe4 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/eos_severe_asthma4","severe_asthma",3,200, "severe_asthma_PB4","severe_asthma","Exp1","severe_asthma")

MNC_ctrl1 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/MNC_healthy1","healthy",3,200, "healthy_MNC1","healhty_MNC","Exp1","healthy")
MNC_ctrl2 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/MNC_healthy2","healthy",3,200, "healthy_MNC2","healhty_MNC","Exp1","healthy")
MNC_ctrl3 <- create_seurat_10X_structured_data( "/scratch/khandl/5.GSE276583/MNC_healthy3","healthy",3,200, "healthy_MNC3","healhty_MNC","Exp1","healthy")

### Merge samples
obj <- merge(PB_ctrl1, y = c(PB_ctrl2, PB_ctrl3,PB_ctrl4,PB_asthma_mild1,PB_asthma_mild2,PB_asthma_mild3,PB_asthma_mild4,
                             PB_asthma_severe1,PB_asthma_severe2,PB_asthma_severe3,PB_asthma_severe4,MNC_ctrl1,MNC_ctrl2,
                             MNC_ctrl3),
             add.cell.ids = c("PB_ctrl1", "PB_ctrl2","PB_ctrl3","PB_ctrl4", "PB_asthma_mild1","PB_asthma_mild2","PB_asthma_mild3",
                              "PB_asthma_mild4","PB_asthma_severe1","PB_asthma_severe2","PB_asthma_severe3","PB_asthma_severe4",
                              "MNC_ctrl1","MNC_ctrl2","MNC_ctrl3"))
obj <- JoinLayers(obj)

### Add mitochondrial percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$tissue)
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

### Rename clusters
current.cluster.ids <- c(0:23)
new.cluster.ids <- c("Eosinophils","lowQ","T","Eosinophils", "lowQ","Eosinophils","?","?",
                     "T","?","Monocytes","PCs","ProNeutro","lowQ","Basophils","B","Mast","?","DCs","lowQ",
                     "Neutrophils","?","?","?")
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
obj$technology <- "10X"
obj$species <- "Hs"
obj$cell_enrichment <- "Eosinophils"
obj$tissue <- "blood"

##### Save object 
saveRDS(obj, file = "/scratch/khandl/technical/seurat_objects/Hs_blood_GSE276583_anno.rds")
