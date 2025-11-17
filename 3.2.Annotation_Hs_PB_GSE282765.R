######### This code integrates and annotates blood data ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_exonic_and_intronic_reads_Hs_PB_healthy_and_CRC.rds")

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

DotPlot(obj, features = unique(c("TPSAB1","KIT","TPSB2", #Mast cells 
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
                                 "COL1A2","PECAM1" ,# Fibroblasts/endothelial 
                                 "CD34","CEBPA","MSI1", #GMPs
                                 "EPX","EAR2","EAR1", #EoP
                                 "ELANE","CEBPE","S100A8" ,#ProNeutro
                                 "LY6C2","CSF1R","IRF8" #ProMono
)) )+ theme(axis.text.x = element_text(angle = 90)) 

### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### nFeature and percent.mito per cluster 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

### rename
current.cluster.ids <- c(0:14)
new.cluster.ids <- c("Eosinophils","Eosinophils","lowQ","lowQ", "ProNeutro","Neutrophils","Neutrophils","T", "Mast_Baso","Mono_Mac","PCs",
                     "Fibroblasts", "Mixed","DCs","?")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn" )

##### subcluster Mono/mac
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Mono_Mac",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mono_Mac_0","Mono_Mac_1"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features = c( "CD300E","EREG","VCAN", #Monocytes
                                    "F13A1", "FOLR2","C1QC" #Macrophages 
), scale = FALSE) + theme(axis.text.x = element_text(angle = 90)) 
# 1 = Macrophages, 0  = Monocytes

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

# rename
current.cluster.ids <- c("?", "DCs", "Eosinophils","Fibroblasts","lowQ","Mast_Baso", "Mixed","Mono_Mac_0","Mono_Mac_1","Neutrophils","PCs","ProNeutro","T")
new.cluster.ids <- c("?", "DCs", "Eosinophils","Fibroblasts","lowQ","Mast_Baso", "Mixed","Monocytes","Macrophages","Neutrophils","PCs","ProNeutro","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster Mast_Baso
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mast_Baso",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mast_Baso_0","Mast_Baso_1","Mast_Baso_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features = unique(c("TPSAB1","KIT","TPSB2", #Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "ICOS","CD8A", "CD3E","TRAC" # T cells 
)), scale = FALSE )+ theme(axis.text.x = element_text(angle = 90)) 

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 2 = Mast cells, 0 = Basophils, 1 = lowQ
VlnPlot(sub_celltype, features = "nFeature_RNA")

# rename
current.cluster.ids <- c( "?", "DCs", "Eosinophils","Fibroblasts","lowQ","Macrophages", 
                          "Mast_Baso_0","Mast_Baso_1","Mast_Baso_2","Mixed","Monocytes","Neutrophils","PCs","ProNeutro","T")
new.cluster.ids <- c( "?", "DCs", "Eosinophils","Fibroblasts","lowQ","Macrophages", 
                      "Basophils","lowQ","Mast","Mixed","Monocytes","Neutrophils","PCs","ProNeutro","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### check annotation
obj <- subCl
Idents(obj) <- "annotation"

DotPlot(obj, features = unique(c("TPSAB1","KIT","TPSB2", #Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,#DCs
                                 "CD300E","EREG","VCAN", #Monocytes
                                 "F13A1", "FOLR2","C1QC" ,#Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", #Neutrophils 
                                 "EPCAM","TFF3","MUC2", #Epithelial 
                                 "COL1A2","PECAM1" ,# Fibroblasts/endothelial 
                                 "CD34","CEBPA","MSI1", #GMPs
                                 "EPX","EAR2","EAR1", #EoP
                                 "ELANE","CEBPE","S100A8" ,#ProNeutro
                                 "LY6C2","CSF1R","IRF8" #ProMono
)) )+ theme(axis.text.x = element_text(angle = 90)) 

obj$sample <- obj$condition
obj$tissue2 <- obj$tissue
obj$tissue <- "blood"
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
