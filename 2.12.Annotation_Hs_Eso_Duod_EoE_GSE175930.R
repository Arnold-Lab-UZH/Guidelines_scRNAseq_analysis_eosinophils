########## This code does sample integration, pre-processing, clustering and annotation of Hs esophagus and duodenum from active EoE and EoE in remission from GSE175930  ##########

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir,"1.2.Functions_Seurat_integration.R"))

##### Seurat object generation 
Duodenum_EoE_rem1 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P249_Duodenum"),"EoE",3,200, "Duodenum_EoE_rem1","Duodenum_EoE_rem","EoE_Exp1","EoE_rem")
Duodenum_EoE_rem2 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P355_Duodenum"),"EoE",3,200, "Duodenum_EoE_rem2","Duodenum_EoE_rem","EoE_Exp2","EoE_rem")
Duodenum_EoE_rem3 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P468_Duodenum"),"EoE",3,200, "Duodenum_EoE_rem3","Duodenum_EoE_rem","EoE_Exp3","EoE_rem")
Duodenum_EoE_rem4 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P392_Duodenum"),"EoE",3,200, "Duodenum_EoE_rem4","Duodenum_EoE_rem","EoE_Exp4","EoE_rem")

Esophagus_EoE_rem1 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P249_Esophagus"),"EoE",3,200, "Esophagus_EoE_rem1","Esophagus_EoE_rem","EoE_Exp1","EoE_rem")
Esophagus_EoE_rem2 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P355_Esophagus"),"EoE",3,200, "Esophagus_EoE_rem2","Esophagus_EoE_rem","EoE_Exp2","EoE_rem")
Esophagus_EoE_rem3 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P468_Esophagus"),"EoE",3,200, "Esophagus_EoE_rem3","Esophagus_EoE_rem","EoE_Exp3","EoE_rem")
Esophagus_EoE_rem4 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P392_Esophagus"),"EoE",3,200, "Esophagus_EoE_rem4","Esophagus_EoE_rem","EoE_Exp4","EoE_rem")

Duodenum_EoE_active1 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P292_Duodenum"),"EoE",3,200, "Duodenum_EoE_active1","Duodenum_EoE_active","EoE_Exp5","EoE_act")
Duodenum_EoE_active2 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P354_Duodenum"),"EoE",3,200, "Duodenum_EoE_active2","Duodenum_EoE_active","EoE_Exp6","EoE_act")
Duodenum_EoE_active3 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P451_Duodenum"),"EoE",3,200, "Duodenum_EoE_active3","Duodenum_EoE_active","EoE_Exp7","EoE_act")
Duodenum_EoE_active5 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P474_Duodenum"),"EoE",3,200, "Duodenum_EoE_active3","Duodenum_EoE_active","EoE_Exp9","EoE_act")
Duodenum_EoE_active6 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P475_Duodenum"),"EoE",3,200, "Duodenum_EoE_active4","Duodenum_EoE_active","EoE_Exp10","EoE_act")

Esophagus_EoE_active1 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P292_Esophagus"),"EoE",3,200, "Esophagus_EoE_active1","Esophagus_EoE_active","EoE_Exp5","EoE_act")
Esophagus_EoE_active2 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P354_Esophagus"),"EoE",3,200, "Esophagus_EoE_active2","Esophagus_EoE_active","EoE_Exp6","EoE_act")
Esophagus_EoE_active3 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P451_Esophagus"),"EoE",3,200, "Esophagus_EoE_active3","Esophagus_EoE_active","EoE_Exp7","EoE_act")
Esophagus_EoE_active4 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P471_Esophagus"),"EoE",3,200, "Esophagus_EoE_active4","Esophagus_EoE_active","EoE_Exp8","EoE_act")
Esophagus_EoE_active5 <- create_seurat_10X_structured_data( file.path(raw_data_GSE175930_dir,"P474_Esophagus"),"EoE",3,200, "Esophagus_EoE_active3","Esophagus_EoE_active","EoE_Exp9","EoE_act")

### Merge samples
tissue <- merge(Duodenum_EoE_rem1, y = c(Duodenum_EoE_rem2,Duodenum_EoE_rem3,Duodenum_EoE_rem4,Esophagus_EoE_rem1,Esophagus_EoE_rem2,Esophagus_EoE_rem3,
                                         Esophagus_EoE_rem4,Duodenum_EoE_active1,Duodenum_EoE_active2,Duodenum_EoE_active3,Duodenum_EoE_active5,Duodenum_EoE_active6,
                                         Esophagus_EoE_active1,Esophagus_EoE_active2,Esophagus_EoE_active3,Esophagus_EoE_active4,Esophagus_EoE_active5),
                add.cell.ids = c("Duodenum_EoE_rem1", "Duodenum_EoE_rem2","Duodenum_EoE_rem3","Duodenum_EoE_rem4", "Esophagus_EoE_rem1","Esophagus_EoE_rem2",
                                 "Esophagus_EoE_rem3","Esophagus_EoE_rem4","Duodenum_EoE_active1","Duodenum_EoE_active2","Duodenum_EoE_active3",
                                 "Duodenum_EoE_active5","Duodenum_EoE_active6","Esophagus_EoE_active1","Esophagus_EoE_active2",
                                 "Esophagus_EoE_active3","Esophagus_EoE_active4","Esophagus_EoE_active5"))
tissue <- JoinLayers(tissue)

### Add mitochondrial percentage per cell 
tissue$percent.mt <- PercentageFeatureSet(tissue, pattern = "^MT-")

### Apply mitochondrial cutoff 
tissue <- subset(tissue, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- tissue
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,  margin = 1,assay = "RNA")
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
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
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
current.cluster.ids <- c(0:14)
new.cluster.ids <- c("T","T","Mast","PCs", "lowQ","Macrophages","Mast","Epithelial","?","?","Eosinophils","Fibroblasts","?","B","?")
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
obj$technology <- "Seq_Well"
obj$species <- "Hs"
obj$cell_enrichment <- "CD45"

current.cluster.ids <- c("Duodenum_EoE_active",     "Duodenum_EoE_rem", "Esophagus_EoE_active",    "Esophagus_EoE_rem")
new.cluster.ids <- c("duodenum",     "duodenum", "esophagus",    "esophagus")
obj$tissue <- plyr::mapvalues(x = obj$tissue, from = current.cluster.ids, to = new.cluster.ids)

##### Save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Hs_esophagus_duodenum_EoE_GSE175930_anno.rds"))

