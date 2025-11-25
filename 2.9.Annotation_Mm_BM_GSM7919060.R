########## This code does sample integration, pre-processing, clustering and annotation of Mm healthy BM from GSM7819060  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Seurat object generation 
BM_eos <- create_seurat_10X_structured_data( "/scratch/khandl/6.GSM7919060/","BM_IL33",3,200, "BM_IL33","BM_IL33","BM_IL33","healthy")

### Add mitochondrial percentage per cell 
BM_eos$percent.mt <- PercentageFeatureSet(BM_eos, pattern = "^mt-")

### Apply mitochondrial cutoff 
BM_eos <- subset(BM_eos, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- BM_eos
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:15, reduction.name = "umap")
DimPlot(obj,reduction = "umap",group.by = "seurat_clusters", label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap",group.by = "condition")
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### nFeature and percent.mito per cluster to exclude low quality clusters 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

### Marker gene expression 
DotPlot(obj, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #B asophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,# DCs
                                  "Cd300e","Ereg","Vcan", # Monocytes
                                  "F13a1", "Folr2","C1qc" ,# Macrophages 
                                  "Spp1","Mmp12","Fn1", # TAMs
                                  "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                                  "Icos","Cd8a", "Cd3e","Trac", # T cells 
                                  "Cd19", "Cr2","Ms4a1", # B cells 
                                  "Jchain","Igha","Igkc", # PCs
                                  "S100a9","S100a8", # Neutrophils 
                                  "Epcam","Tff3","Muc2", # Epithelial 
                                  "Col1a1", # Fibroblasts/
                                  "Pecam1" ,# Endothelial 
                                  "Mki67","Msi1","Meis1", # MPPs
                                  "Cd34","Cebpa", # GMPs
                                  "Elane","Cepbe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:6)
new.cluster.ids <- c("lowQ","ProNeutro","EoP","GMPs", "lowQ","Mixed","Neutrophils")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,reduction = "umap" )

##### Check annotation
Idents(obj) <- "annotation"
DotPlot(obj, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #B asophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,# DCs
                                  "Cd300e","Ereg","Vcan", # Monocytes
                                  "F13a1", "Folr2","C1qc" ,# Macrophages 
                                  "Spp1","Mmp12","Fn1", # TAMs
                                  "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                                  "Icos","Cd8a", "Cd3e","Trac", # T cells 
                                  "Cd19", "Cr2","Ms4a1", # B cells 
                                  "Jchain","Igha","Igkc", # PCs
                                  "S100a9","S100a8", # Neutrophils 
                                  "Epcam","Tff3","Muc2", # Epithelial 
                                  "Col1a1", # Fibroblasts/
                                  "Pecam1" ,# Endothelial 
                                  "Mki67","Msi1","Meis1", # MPPs
                                  "Cd34","Cebpa", # GMPs
                                  "Elane","Cepbe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
obj$technology <- "10X"
obj$species <- "Mm"
obj$cell_enrichment <- "Eosinophils"
obj$tissue <- "bm"

##### Save object 
saveRDS(obj, file = "/scratch/khandl/technical/seurat_objects/Mm_bm_GSM7919060_anno.rds")
