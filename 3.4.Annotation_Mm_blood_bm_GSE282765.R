########## This code integrates colon data from mouse data blood and bm tumor and non-tumor ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds")

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
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label= TRUE, label.size = 8)
obj <- JoinLayers(obj)

### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### nFeature and percent.mito per cluster 
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
VlnPlot(obj, features = "percent.mt", pt.size = 0)

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

# rename
current.cluster.ids <- c(0:22)
new.cluster.ids <- c("Neutrophils","Neutrophils","Neutrophils","Neutrophils", "Monocytes", "lowQ",
                     "?","ProMono_GMPs","ProNeutro","?","ProMono_GMPs",
                     "T_PCs","Eosinophils","?","?","Basophils",
                     "Monocytes","Mixed", "?","DCs","EoP",
                     "?","Macrophages")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster ProMono_GMPs
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "ProMono_GMPs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "ProMono_GMPs_0","ProMono_GMPs_1","ProMono_GMPs_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #Basophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,#DCs
                                  "Cd300e","Ereg","Vcan", #Monocytes
                                  "F13a1", "Folr2","C1qc" ,#Macrophages 
                                  "Mki67","Msi1","Meis1", #MPPs
                                  "Cd34","Cebpa", #GMPs
                                  "Elane","Cepbe", #ProNeutro
                                  "Ly6c2","Csf2r","Irf8", #ProMono
                                  "Flt3","Dntt","Il7r", #CLPs
                                  "Epx","Ear1","Ear2", #EoP
                                  "Neat1","Xist","Malat1" #cytoplasma lost 
)), scale = FALSE) + theme(axis.text.x = element_text(angle = 90)) 
# 1 = GMPs, 0,2 = ProMono

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# rename
current.cluster.ids <- c( "?","Basophils", "DCs","EoP","Eosinophils","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                          "ProMono_GMPs_0","ProMono_GMPs_1","ProMono_GMPs_2", "ProNeutro","T_PCs")
new.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                     "ProMono","GMPs","ProMono", "ProNeutro","T_PCs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### subcluster T_PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T_PCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_PCs_0","T_PCs_1","T_PCs_2","T_PCs_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 1 = T ; 2 PCs 
FeaturePlot(sub_celltype, features = c("Cd19","Cd3e","Mpl","Mki67","Meis1"))
# Mpl is a marker for hemaptopoietic stem cells 

# rename
current.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","GMPs","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                         "ProMono","ProNeutro", "T_PCs_0","T_PCs_1","T_PCs_2","T_PCs_3")
new.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","GMPs","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                     "ProMono","ProNeutro", "Mixed","T","B","HSCs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### check annotation
obj <- subCl
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
                                  "Col1a2","Pecam1" ,# Fibroblasts/endothelial 
                                  "Mpl","Meis1", #HPCs
                                  "Cd34","Cebpa", #GMPs
                                  "Elane","Cepbe", #ProNeutro
                                  "Ly6c2","Csf2r","Irf8", #ProMono
                                  "Flt3","Dntt","Il7r", #CLPs
                                  "Epx","Ear1","Ear2", #EoP
                                  "Neat1","Xist","Malat1" #cytoplasma lost 
))) + theme(axis.text.x = element_text(angle = 90)) 

obj$sample <- obj$condition

current.cluster.ids <- c("blood_tumor","blood_wt","bm_tumor","bm_wt")
new.cluster.ids <- c("blood","blood","bm","bm")
obj$tissue <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")

