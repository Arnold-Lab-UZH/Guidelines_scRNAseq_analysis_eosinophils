########## This code does sample integration, pre-processing, clustering and annotation of Mm Il5tg tissues from GGSE182001 ##########

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir,"1.2.Functions_Seurat_integration.R"))

##### Seurat object generation 
### forced cell determination - intronic and exonic reads 
stomach <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_stomach_ST01_Expression_Data.st"), 
  project = "steady_state", condition = "stomach",3,200)

colon <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_colon_ST02_Expression_Data.st"), 
  project = "steady_state", condition = "colon",3,200)

small_int <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_SI_ST03_Expression_Data.st"), 
  project = "steady_state", condition = "small_int",3,200)

spleen <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_spleen_ST04_Expression_Data.st"), 
  project = "steady_state", condition = "spleen",3,200)

blood <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_blood_ST08_Expression_Data.st"), 
  project = "steady_state", condition = "blood",3,200)

bm <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE182001_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_bm_ST06_Expression_Data.st"), 
  project = "steady_state", condition = "bm",3,200)

### Merge samples
merged <- merge(stomach, y = c(colon,small_int,spleen,  blood, bm),
                add.cell.ids = c("stomach","colon","SI","spleen", "blood","bm"))
merged <- JoinLayers(merged)

### Add mitochondrial percentage per cell 
merged$percent.mt <- PercentageFeatureSet(merged, pattern = "^mt-")

### Add conditions to metadata 
merged$cell_determination <- "forced"
merged$reads <- "intronic_and_exonic"
merged$species <- "Mm"
merged$technology <- "BD_Rhapsody"
merged$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(merged, file = file.path(seurat_objects_dir,"Forced_cell_determination_intronic_and_exonic_reads_Mm_il5tg_steady_state.rds"))

###### load R object 
obj <- readRDS(file = file.path(seurat_objects_dir,"Forced_cell_determination_intronic_and_exonic_reads_Mm_il5tg_steady_state.rds"))

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,  margin = 1,assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15)
DimPlot(obj,reduction = "umap",group.by = "seurat_clusters", label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
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
                                  "Elane","Cebpe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:15)
new.cluster.ids <- c("Eosinophils","Eosinophils","Eosinophils","Macrophages", "EoP", "lowQ",
                     "?","Mix_progenitor", "Endothelial_Fibroblasts","T_PCs","Neutrophils",
                     "?","ProNeutro","Mix_progenitor","Epithelial","?")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE)

### Subcluster Mix_progenitor
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Mix_progenitor",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mix_progenitor_0","Mix_progenitor_1","Mix_progenitor_2","Mix_progenitor_3"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

DotPlot(sub_celltype, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
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
                                  "Elane","Cebpe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 
# 3 = PCs, everything else mixed 

# Rename
current.cluster.ids <- c( "?","Endothelial_Fibroblasts","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                          "Mix_progenitor_0","Mix_progenitor_1","Mix_progenitor_2","Mix_progenitor_3",
                          "Neutrophils","ProNeutro", "T_PCs")
new.cluster.ids <- c( "?","Endothelial_Fibroblasts","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                      "Mixed","Mixed","Mixed","PCs",
                      "Neutrophils","ProNeutro", "T_PCs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### Subcluster T_PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T_PCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_PCs_0","T_PCs_1"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

DotPlot(sub_celltype, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
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
                                           "Elane","Cebpe", # ProNeutro
                                           "Ly6c2","Csf2r","Irf8", # ProMono
                                           "Flt3","Dntt","Il7r", # CLPs
                                           "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 
# everything else mixed 

# Rename
current.cluster.ids <- c(  "?","Endothelial_Fibroblasts","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                           "Mixed", "Neutrophils","PCs", "ProNeutro", "T_PCs_0","T_PCs_1")
new.cluster.ids <- c( "?","Endothelial_Fibroblasts","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                      "Mixed", "Neutrophils","PCs", "ProNeutro", "Mixed","Mixed")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

### Subcluster Fib_Endo
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Endothelial_Fibroblasts",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.05)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Endothelial_Fibroblasts_0","Endothelial_Fibroblasts_1"))
DimPlot(sub_celltype, label = TRUE, group.by = "sub.cluster")

FeaturePlot(sub_celltype, features = c("Pecam1","Col1a2"))
# 0 = Endothelial, 1 = Fibroblasts

# Rename
current.cluster.ids <- c(  "?","Endothelial_Fibroblasts_0","Endothelial_Fibroblasts_1","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                           "Mixed", "Neutrophils","PCs", "ProNeutro")
new.cluster.ids <- c( "?","Endothelial","Fibroblasts","EoP","Eosinophils","Epithelial","lowQ", "Macrophages",
                      "Mixed", "Neutrophils","PCs", "ProNeutro")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE)

##### Check annotation
obj <- subCl
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
                                  "Elane","Cebpe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
))) + theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
obj$sample <- obj$condition
obj$tissue <- obj$condition
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj, file.path(seurat_objects_dir,"Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds"))
