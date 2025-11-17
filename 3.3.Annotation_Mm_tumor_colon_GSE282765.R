########## This code integrates colon data from mouse data tumor and non-tumor ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_intronic_and_exonic_reads_Mm_CRC_AKPS_tumor_NAT_diss.rds")

obj <- subset(obj, subset = percent.mt < 25)

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

#### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)                     
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
obj <- JoinLayers(obj)

### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### nFeature and percent.mito per cluster 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

DotPlot(obj, features =  c("Tpsab1", "Tpsb2","Kit",#Mast cells 
                           "Mcpt8","Cd200r3","Clec12a", #Basophils
                           "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                           "Clec9a","Flt3","Xcr1" ,#DCs
                           "Cd300e","Ereg","Vcan", #Monocytes
                           "F13a1", "Folr2","C1qc" ,#Macrophages 
                           "Spp1","Mmp12","Fn1", # TAMs
                           "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                           "Icos","Cd8a", "Cd3e","Trac", # T cells 
                           "Cd19", "Cr2","Ms4a1", # B cells 
                           "Igha","Jchain","Igkc", # PCs
                           "S100a8","S100a9", #Neutrophils 
                           "Epcam","Tff3","Muc2", #Epithelial 
                           "Col1a2","Pecam1" # Fibroblasts/endothelial 
)) + theme(axis.text.x = element_text(angle = 90)) 

# rename
current.cluster.ids <- c(0:21)
new.cluster.ids <- c("Neutrophils","T","TAMs","lowQ", "B", "Neutrophils",
                     "Macrophages","Eosinophils","PCs","Baso_Mast_mixed","lowQ",
                     "?","DCs","Mixed","PCs","Mixed",
                     "Mixed","Mixed","Mixed","?","?","Mixed")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster Baso_Mast_mixed
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

# rename
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

##### subcluster TAMs to find monoytes 
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

# rename
current.cluster.ids <- c(  "?","B", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Neutrophils","PCs","T", "TAMs_0","TAMs_1")
new.cluster.ids <- c( "?","B", "DCs", "Eosinophils","lowQ",
                      "Macrophages","Mast", "Mixed","Neutrophils","PCs","T", "TAMs","Monocytes")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster T 
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

# rename
current.cluster.ids <- c(  "?","B", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes", "Neutrophils","PCs","T_0","T_1","T_2","T_3","TAMs")
new.cluster.ids <- c(   "?","B", "DCs", "Eosinophils","lowQ",
                        "Macrophages","Mast", "Mixed","Monocytes", "Neutrophils","PCs","T","T","T","Mixed","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster B
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

# rename
current.cluster.ids <- c(  "?", "B_0","B_1","B_2", "DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
new.cluster.ids <- c(   "?", "B","B","lowQ", "DCs", "Eosinophils","lowQ",
                        "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster PCs
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

# rename
current.cluster.ids <- c(  "?", "B","DCs", "Eosinophils","lowQ",
                           "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs_0","PCs_1","PCs_2","PCs_3","PCs_4","PCs_5","T","TAMs")
new.cluster.ids <- c(  "?", "B","DCs", "Eosinophils","lowQ",
                       "Macrophages","Mast", "Mixed","Monocytes","Neutrophils","PCs","PCs","PCs","Mixed","PCs","Mixed","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### check annotation
obj <- subCl
Idents(obj) <- "annotation"
DotPlot(obj, features =  c("Tpsab1", "Tpsb2","Kit", # Mast cells
                           "Mcpt8","Cd200r3","Clec12a", #Basophils
                           "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                           "Clec9a","Flt3","Xcr1" ,#DCs
                           "Cd300e","Ereg","Vcan", #Monocytes
                           "F13a1", "Folr2","C1qc" ,#Macrophages 
                           "Spp1","Mmp12","Fn1", # TAMs
                           "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                           "Icos","Cd8a", "Cd3e","Trac", # T cells 
                           "Cd19", "Cr2","Ms4a1", # B cells 
                           "Jchain","Igkc", # PCs
                           "S100a8","S100a9", #Neutrophils 
                           "Epcam","Tff3","Muc2", #Epithelial 
                           "Col1a2","Pecam1" # Fibroblasts/endothelial 
)) + theme(axis.text.x = element_text(angle = 90)) 

obj$sample <- obj$condition

current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids <- c("colon","colon","tumor","tumor")
obj$tissue <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
