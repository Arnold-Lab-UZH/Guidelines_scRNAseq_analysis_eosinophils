##########  Pre-processing and annotation of mouse BM Siglecf sorted eosinophils GSM7919060 ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

BM_eos <- create_seurat_public_data( "/scratch/khandl/6.GSM7919060/","BM_IL33",3,200, "BM_IL33","BM_IL33","BM_IL33","healthy")

### Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
BM_eos$percent.mt <- PercentageFeatureSet(BM_eos, pattern = "^mt-")

### apply mitochondrial and nFeature cutoffs
BM_eos <- subset(BM_eos, subset = percent.mt < 25)

##### Annotation 
obj <- BM_eos
### pre-processing
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

DotPlot(obj, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #Basophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,#DCs
                                  "Cd300e","Ereg","Vcan", #Monocytes
                                  "F13a1", "Folr2","C1qc" ,#Macrophages 
                                  "Spp1","Mmp12","Fn1", # TAMs
                                  "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                                  "Icos","Cd8a", "Cd3e","Trac", # T cells 
                                  "Cd19", "Cr2","Ms4a1", # B cells 
                                  "Jchain","Igha","Igkc", # PCs
                                  "S100a9","S100a8", #Neutrophils 
                                  "Epcam","Tff3","Muc2", #Epithelial 
                                  "Col1a1","Pecam1" ,# Fibroblasts/endothelial 
                                  "Mki67","Msi1","Meis1", #MPPs
                                  "Cd34","Cebpa", #GMPs
                                  "Elane","Cepbe", #ProNeutro
                                  "Ly6c2","Csf2r","Irf8", #ProMono
                                  "Flt3","Dntt","Il7r", #CLPs
                                  "Epx","Ear1","Ear2", #EoP
                                  "Neat1","Xist","Malat1" #cytoplasma lost 
))) + theme(axis.text.x = element_text(angle = 90)) 

### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### nFeature and percent.mito per cluster 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

### rename
current.cluster.ids <- c(0:6)
new.cluster.ids <- c("lowQ","ProNeutro","EoP","GMPs", "lowQ","Mixed","Neutrophils")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,reduction = "umap" )

##### check annotataion
Idents(obj) <- "annotation"
DotPlot(obj, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #Basophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,#DCs
                                  "Cd300e","Ereg","Vcan", #Monocytes
                                  "F13a1", "Folr2","C1qc" ,#Macrophages 
                                  "Spp1","Mmp12","Fn1", # TAMs
                                  "Siglecf","Alox15","Ccr3","F5", # Eosinophils 
                                  "Icos","Cd8a", "Cd3e","Trac", # T cells 
                                  "Cd19", "Cr2","Ms4a1", # B cells 
                                  "Jchain","Igha","Igkc", # PCs
                                  "S100a9","S100a8", #Neutrophils 
                                  "Epcam","Tff3","Muc2", #Epithelial 
                                  "Col1a1","Pecam1" ,# Fibroblasts/endothelial 
                                  "Mki67","Msi1","Meis1", #MPPs
                                  "Cd34","Cebpa", #GMPs
                                  "Elane","Cepbe", #ProNeutro
                                  "Ly6c2","Csf2r","Irf8", #ProMono
                                  "Flt3","Dntt","Il7r", #CLPs
                                  "Epx","Ear1","Ear2", #EoP
                                  "Neat1","Xist","Malat1" #cytoplasma lost 
))) + theme(axis.text.x = element_text(angle = 90)) 

obj$technology <- "10X"
obj$species <- "Mm"
obj$cell_enrichment <- "Eosinophils"
obj$tissue <- "bm"

### save object 
saveRDS(obj, file = "/scratch/khandl/technical/seurat_objects/Mm_bm_GSM7919060_anno.rds")
