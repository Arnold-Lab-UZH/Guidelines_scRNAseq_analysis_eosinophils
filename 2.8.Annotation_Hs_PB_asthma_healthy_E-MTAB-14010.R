########## This code does sample integration, pre-processing, clustering and annotation of Hs PB from asthma patients and healthy individuals from E-MTAB-14010  ##########

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir,"1.2.Functions_Seurat_integration.R"))

##### Seurat object generation 
PB_ctrl1 <- create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Control1"),"Control1",3,200, "PB_control1","PB_healthy","Exp1","healthy")
PB_ctrl2 <- create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Control2"),"Control",3,200, "PB_control2","PB_healthy","Exp1","healthy")
PB_ctrl3 <- create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Control3"),"Control",3,200, "PB_control3","PB_healthy","Exp1","healthy")

PB_asthma1 <- create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Asthma1"),"Asthma",3,200, "PB_asthma1","PB_Asthma","Exp1","asthma")
PB_asthma2 <-  create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Asthma2"),"Asthma",3,200, "PB_asthma2","PB_Asthma","Exp1","asthma")
PB_asthma3 <- create_seurat_10X_structured_data( file.path(raw_data_E_MTAB_14010_dir,"Asthma3"),"Asthma",3,200, "PB_asthma3","PB_Asthma","Exp1","asthma")

### Merge samples
obj <- merge(PB_ctrl1, y = c(PB_ctrl2, PB_ctrl3,PB_asthma1,PB_asthma2,PB_asthma3),
               add.cell.ids = c("PB_ctrl1", "PB_ctrl2","PB_ctrl3","PB_asthma1", "PB_asthma2","PB_asthma3"))
obj <- JoinLayers(obj)

### Add mitochondrial percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.1, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:20, reduction.name = "umap")
DimPlot(obj,reduction = "umap",group.by = "seurat_clusters", label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap",group.by = "condition")
obj <- JoinLayers(obj)

##### Cluster annotation 
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
# 4 = T, 3 = ProNeutro, 0,1,2 = Eosinophils, 5 = Neutrophils
FeaturePlot(obj, features = "CD3E")

### Rename clusters
current.cluster.ids <- c(0:5)
new.cluster.ids <- c("Eosinophils","Eosinophils","Eosinophils","ProNeutro","Mixed","Neutrophils")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,reduction = "umap" )

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
obj$technology <- "BD"
obj$species <- "Hs"
obj$cell_enrichment <- "Eosinophils"
obj$tissue <- "blood"

##### Save object 
saveRDS(obj, file = file.path(seurat_objects_dir,"Hs_blood_E-MTAB-14010_anno.rds"))

