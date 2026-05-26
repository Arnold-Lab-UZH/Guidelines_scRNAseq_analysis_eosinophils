########## This code does sample integration, pre-processing, clustering and annotation of Mm eosinophil subtypes from GSE182001  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load annotated R object from 2.5 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Extract Eosinophils and EoP 
Idents(obj) <- "annotation"
obj <- subset(obj, idents =c("EoP","Eosinophils"))

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
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:20, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 10)
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### nFeature and percent.mito per cluster to exclude low quality clusters 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")

### Marker gene expression 
FeaturePlot(obj, features = c("Siglecf","Ccr3","Syne1","Alox15"),reduction = "umap.mnn")
FeaturePlot(obj, features = c("Epx","Ear1","Alox15","Aldh2","Ccl9","Retnla","Icosl","Il4","Cd80","Cd274"),reduction = "umap.mnn")

goi <- c("Mki67","Epx","Prg3", #progenitor
         "Ear1","Ear2","Alox15","Cebpe", # immature 
         "Aldh2","S100a6","Retnla","Ccl9","Il1rl1","Cd24a", # circulating 
         "Mmp9","Icosl","Il4","Tgfb1","Pirb","Rara", # basal 
         "Cd80","Cd274","Ptgs2","Il1rn","Il1b","Vegfa" # active 
         )
Idents(obj) <- "mnn.clusters"
DotPlot(obj, features = goi,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

### Rename
current.cluster.ids <- c(0:6)
new.cluster.ids <- c("circulating","lowQ","active","basal","basal" ,"immature","progenitor")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Check annotation 
Idents(obj) <- "annotation"
DotPlot(obj, features = goi,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

### Remove lowQ 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("circulating","active","basal", "progenitor", "immature"))

##### Add conditions to metadata 
obj$sample <- obj$condition
obj$tissue <- obj$condition
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj,"/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_annotation.rds")

