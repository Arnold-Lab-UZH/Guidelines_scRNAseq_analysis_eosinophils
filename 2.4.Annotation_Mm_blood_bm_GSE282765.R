########## This code does sample integration, pre-processing, clustering and annotation of blood and BM of healthy mice and AKPS CRC mice from GSE282765  ##########

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir,"1.2.Functions_Seurat_integration.R"))

##### Seurat object generation 
### Forced cell determination - intronic and exonic reads 
blood_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_blood_healthy_ST11_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

blood_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_blood_AKPS_tumor_ST12_Expression_Data.st"), 
  project = "blood_tumor", condition = "blood_tumor",3,200)

bm_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_bm_healthy_ST09_Expression_Data.st"), 
  project = "bm_wt", condition = "bm_wt",3,200)

bm_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_forced_intron_exon_dir, "Forced_cell_determination_intronic_and_exonic_Mm_bm_AKPS_tumor_ST10_Expression_Data.st"), 
  project = "bm_tumor", condition = "bm_tumor",3,200)

### Merge samples
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt-")

### Add conditions to metadata 
tumor$cell_determination <- "forced"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(tumor, file = file.path(seurat_objects_dir,"Forced_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds"))

##### Load R object 
obj <- readRDS(file = file.path(seurat_objects_dir,"Forced_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds"))

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)                     
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label= TRUE, label.size = 8)
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
current.cluster.ids <- c(0:22)
new.cluster.ids <- c("Neutrophils","Neutrophils","Neutrophils","Neutrophils", "Monocytes", "lowQ",
                     "?","ProMono_GMPs","ProNeutro","?","ProMono_GMPs",
                     "T_PCs","Eosinophils","?","?","Basophils",
                     "Monocytes","Mixed", "?","DCs","EoP",
                     "?","Macrophages")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster ProMono_GMPs
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "ProMono_GMPs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "ProMono_GMPs_0","ProMono_GMPs_1","ProMono_GMPs_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features =  unique(c("Tpsab1", "Tpsb2","Kit", # Mast cells
                                  "Mcpt8","Cd200r3","Clec12a", #B asophils
                                  "Cpa3","Fcer1a","Ms4a2","Hdc", #Mast and Basophils
                                  "Clec9a","Flt3","Xcr1" ,# DCs
                                  "Cd300e","Ereg","Vcan", # Monocytes
                                  "F13a1", "Folr2","C1qc" ,# Macrophages 
                                  "Mki67","Msi1","Meis1", # MPPs
                                  "Cd34","Cebpa", # GMPs
                                  "Elane","Cebpe", # ProNeutro
                                  "Ly6c2","Csf2r","Irf8", # ProMono
                                  "Flt3","Dntt","Il7r", # CLPs
                                  "Epx","Ear1","Ear2" # EoP
)), scale = FALSE) + theme(axis.text.x = element_text(angle = 90)) 
# 1 = GMPs, 0,2 = ProMono

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

# Rename
current.cluster.ids <- c( "?","Basophils", "DCs","EoP","Eosinophils","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                          "ProMono_GMPs_0","ProMono_GMPs_1","ProMono_GMPs_2", "ProNeutro","T_PCs")
new.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                     "ProMono","GMPs","ProMono", "ProNeutro","T_PCs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster T_PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T_PCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_PCs_0","T_PCs_1","T_PCs_2","T_PCs_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 1 = T ; 2 PCs 
FeaturePlot(sub_celltype, features = c("Cd19","Cd3e","Mpl","Mki67","Meis1"))
# Mpl is a marker for hemaptopoietic stem cells 

# Rename
current.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","GMPs","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                         "ProMono","ProNeutro", "T_PCs_0","T_PCs_1","T_PCs_2","T_PCs_3")
new.cluster.ids <- c("?","Basophils", "DCs","EoP","Eosinophils","GMPs","lowQ","Macrophages","Mixed","Monocytes","Neutrophils",
                     "ProMono","ProNeutro", "Mixed","T","B","HSCs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Check annotation
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
current.cluster.ids <- c("blood_tumor","blood_wt","bm_tumor","bm_wt")
new.cluster.ids <- c("blood","blood","bm","bm")
obj$tissue <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj, file.path(seurat_objects_dir,"Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))

