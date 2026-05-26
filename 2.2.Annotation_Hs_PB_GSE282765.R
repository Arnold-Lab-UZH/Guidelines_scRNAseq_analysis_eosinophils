########## This code does sample integration, pre-processing, clustering and annotation of Hs PB healthy and CRC from GSE282765  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Seurat object generation 
### Forced cell determination intronic and exonic reads 
## Healthy individuals 
H1 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H1_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                            "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H2_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                            "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H3_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                            "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H4_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                            "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H5_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                            "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H6_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                            "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")
H7 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H7_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                            "H7",3,200,  "H7_blood","blood_healthy","Exp1","healthy")
H8 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H8_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                            "H8",3,200,  "H8_blood","blood_healthy","Exp2","healthy")
H9 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H9_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                            "H9",3,200,  "H9_blood","blood_healthy","Exp3","healthy")

## CRC patients 
P1_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P1_Hs_PB_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                  "P1",3,200,  "P1_blood","blood_patient","Exp1","patient")
P2_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_PB_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P2",3,200,  "P2_blood","blood_patient","Exp2","patient")
P3_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_PB_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                  "P3",3,200,  "P3_blood","blood_patient","Exp3","patient")
P4_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_PB_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P4",3,200,  "P4_blood","blood_patient","Exp4","patient")
P5_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_PB_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                  "P5",3,200,  "P5_blood","blood_patient","Exp5","patient")
P6_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_PB_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                  "P6",3,200,  "P6_blood","blood_patient","Exp7","patient")
P7_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_PB_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                  "P7",3,200,  "P7_blood","blood_patient","Exp8","patient")

## Merge samples
blood <- merge(H1, y = c(H2,H3,H4,H5,H6,H7,H8,H9,P1_blood,P2_blood,P3_blood,P4_blood,P5_blood,P6_blood,P7_blood),
               add.cell.ids = c("b1", "b2","b3","b4", "b5","b6","b7","b8","b9","pb1", "pb2","pb3","pb4", "pb5","pb6","pb7"))
blood <- JoinLayers(blood)

### Add mitochondrial percentage per cell 
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

### Add conditions to metadata 
blood$cell_determination <- "forced"
blood$reads <- "intronic_and_exonic"
blood$species <- "Hs"
blood$technology <- "BD_Rhapsody"
blood$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(blood, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_and_intronic_reads_Hs_PB_healthy_and_CRC.rds")

##### Load R object 
obj <- readRDS(file = "/scratch/khandl/4.Technical/Forced_cell_determination_exonic_and_intronic_reads_Hs_PB_healthy_and_CRC.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
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
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
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

### Rename clusters
current.cluster.ids <- c(0:14)
new.cluster.ids <- c("Eosinophils","Eosinophils","lowQ","lowQ", "ProNeutro","Neutrophils","Neutrophils","T", "Mast_Baso","Mono_Mac","PCs",
                     "Fibroblasts", "Mixed","DCs","?")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn" )

### Subcluster Mono/mac
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Mono_Mac",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mono_Mac_0","Mono_Mac_1"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features = c( "CD300E","EREG","VCAN", #Monocytes
                                    "F13A1", "FOLR2","C1QC" # Macrophages 
), scale = FALSE) + theme(axis.text.x = element_text(angle = 90)) 
# 1 = Macrophages, 0  = Monocytes

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "DCs", "Eosinophils","Fibroblasts","lowQ","Mast_Baso", "Mixed","Mono_Mac_0","Mono_Mac_1","Neutrophils","PCs","ProNeutro","T")
new.cluster.ids <- c("?", "DCs", "Eosinophils","Fibroblasts","lowQ","Mast_Baso", "Mixed","Monocytes","Macrophages","Neutrophils","PCs","ProNeutro","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Subcluster Mast_Baso
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mast_Baso",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mast_Baso_0","Mast_Baso_1","Mast_Baso_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DotPlot(sub_celltype, features = unique(c("TPSAB1","KIT","TPSB2", #M ast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "ICOS","CD8A", "CD3E","TRAC" # T cells 
)), scale = FALSE )+ theme(axis.text.x = element_text(angle = 90)) 

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
# 2 = Mast cells, 0 = Basophils, 1 = lowQ
VlnPlot(sub_celltype, features = "nFeature_RNA")

# Rename
current.cluster.ids <- c( "?", "DCs", "Eosinophils","Fibroblasts","lowQ","Macrophages", 
                          "Mast_Baso_0","Mast_Baso_1","Mast_Baso_2","Mixed","Monocytes","Neutrophils","PCs","ProNeutro","T")
new.cluster.ids <- c( "?", "DCs", "Eosinophils","Fibroblasts","lowQ","Macrophages", 
                      "Basophils","lowQ","Mast","Mixed","Monocytes","Neutrophils","PCs","ProNeutro","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### Check annotation
obj <- subCl
Idents(obj) <- "annotation"

DotPlot(obj, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
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
obj$sample <- obj$condition
obj$tissue2 <- obj$tissue
obj$tissue <- "blood"
obj$technology <- "BD"
obj$enrichment <- obj$cell_enrichment

##### Save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
