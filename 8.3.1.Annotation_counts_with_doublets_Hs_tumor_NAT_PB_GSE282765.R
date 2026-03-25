########## This code does annotation of count matrices including multiplets before cell hashing deconvolution   ##########
# Hs tumor, NAT, PB data from GSE282765

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Pre-processing: Loading of unfiltered matrix, EmptyDrops and merging of resulting data 
### P1 
## load data 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P1_tumor_NAT_blood_Expression_Data_Unfiltered.st")

## run EmptyDrops
set.seed(100)
# Test ambient means that cells with lower than 100 UMI counts are also tested against the ambient profile 
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
# extract real cells ids 
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)

## extract counts from only real cells 
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]

## generate Seurat object from real cells 
P1_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P1",3,200,  "P1","with_multiplet","Exp1","patient")

## add singlet/multiplet condition
P1 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P1_Sample_Tag_Calls.csv",skip = 7)
P1_multiplets <- P1[P1$Sample_Tag %in% "Multiplet",]$Cell_Index
P1_singlets <- P1[P1$Sample_Tag %in% c("SampleTag01_hs","SampleTag02_hs","SampleTag03_hs","SampleTag04_hs"),]$Cell_Index

P1_seurat$cell_hashing_cond <- NA
P1_seurat@meta.data <- P1_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P1_seurat@meta.data) %in% P1_multiplets ~ "Doublet",
    rownames(P1_seurat@meta.data) %in% P1_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P1_seurat) <- "cell_hashing_cond"
P1_seurat <- subset(P1_seurat, idents = c("Doublet","Singlet"))

### P2 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P2_tumor_NAT_blood_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P2_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P2",3,200,  "P2","with_multiplet","Exp2","patient")

## add singlet/multiplet condition
P2 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P2_Sample_Tag_Calls.csv",skip = 7)
P2_multiplets <- P2[P2$Sample_Tag %in% "Multiplet",]$Cell_Index
P2_singlets <- P2[P2$Sample_Tag %in% c("SampleTag05_hs","SampleTag06_hs","SampleTag07_hs","SampleTag08_hs"),]$Cell_Index

P2_seurat$cell_hashing_cond <- NA
P2_seurat@meta.data <- P2_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P2_seurat@meta.data) %in% P2_multiplets ~ "Doublet",
    rownames(P2_seurat@meta.data) %in% P2_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P2_seurat) <- "cell_hashing_cond"
P2_seurat <- subset(P2_seurat, idents = c("Doublet","Singlet"))

## P3 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P3_tumor_NAT_blood_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P3_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P3",3,200,  "P3","with_multiplet","Exp3","patient")

## add singlet/multiplet condition
P3 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P3_Sample_Tag_Calls.csv",skip = 7)
P3_multiplets <- P3[P3$Sample_Tag %in% "Multiplet",]$Cell_Index
P3_singlets <- P3[P3$Sample_Tag %in% c("SampleTag09_hs","SampleTag10_hs","SampleTag11_hs","SampleTag12_hs"),]$Cell_Index

P3_seurat$cell_hashing_cond <- NA
P3_seurat@meta.data <- P3_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P3_seurat@meta.data) %in% P3_multiplets ~ "Doublet",
    rownames(P3_seurat@meta.data) %in% P3_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P3_seurat) <- "cell_hashing_cond"
P3_seurat <- subset(P3_seurat, idents = c("Doublet","Singlet"))

## P4 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P4_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P4_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P4",3,200,  "P4","with_multiplet","Exp4","patient")

## add singlet/multiplet condition
P4 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P4_Sample_Tag_Calls.csv",skip = 7)
P4_multiplets <- P4[P4$Sample_Tag %in% "Multiplet",]$Cell_Index
P4_singlets <- P4[P4$Sample_Tag %in% c("SampleTag01_hs","SampleTag02_hs","SampleTag03_hs","SampleTag04_hs"),]$Cell_Index

P4_seurat$cell_hashing_cond <- NA
P4_seurat@meta.data <- P4_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P4_seurat@meta.data) %in% P4_multiplets ~ "Doublet",
    rownames(P4_seurat@meta.data) %in% P4_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P4_seurat) <- "cell_hashing_cond"
P4_seurat <- subset(P4_seurat, idents = c("Doublet","Singlet"))

## P5 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P5_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P5_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P5",3,200,  "P5","with_multiplet","Exp5","patient")

## add singlet/multiplet condition
P5 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P5_Sample_Tag_Calls.csv",skip = 7)
P5_multiplets <- P5[P5$Sample_Tag %in% "Multiplet",]$Cell_Index
P5_singlets <- P5[P5$Sample_Tag %in% c("SampleTag06_hs","SampleTag07_hs","SampleTag08_hs"),]$Cell_Index

P5_seurat$cell_hashing_cond <- NA
P5_seurat@meta.data <- P5_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P5_seurat@meta.data) %in% P5_multiplets ~ "Doublet",
    rownames(P5_seurat@meta.data) %in% P5_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P5_seurat) <- "cell_hashing_cond"
P5_seurat <- subset(P5_seurat, idents = c("Doublet","Singlet"))

## P6 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P6_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P6_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P6",3,200,  "P6","with_multiplet","Exp7","patient")

## add singlet/multiplet condition
P6 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P6_Sample_Tag_Calls.csv",skip = 7)
P6_multiplets <- P6[P6$Sample_Tag %in% "Multiplet",]$Cell_Index
P6_singlets <- P6[P6$Sample_Tag %in% c("SampleTag02_hs","SampleTag05_hs","SampleTag06_hs","SampleTag07_hs"),]$Cell_Index

P6_seurat$cell_hashing_cond <- NA
P6_seurat@meta.data <- P6_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P6_seurat@meta.data) %in% P6_multiplets ~ "Doublet",
    rownames(P6_seurat@meta.data) %in% P6_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P6_seurat) <- "cell_hashing_cond"
P6_seurat <- subset(P6_seurat, idents = c("Doublet","Singlet"))

## P7 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P7_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
real_cells_id <- rownames(df)
counts_data <- counts_data[,colnames(counts_data) %in% real_cells_id]
P7_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P7",3,200,  "P7","with_multiplet","Exp8","patient")

## add singlet/multiplet condition
P7 <- read.csv("/scratch/khandl/Technical_count_matrices/Sample_tag_calls/Hs_P7_Sample_Tag_Calls.csv",skip = 7)
P7_multiplets <- P7[P7$Sample_Tag %in% "Multiplet",]$Cell_Index
P7_singlets <- P7[P7$Sample_Tag %in% c("SampleTag10_hs","SampleTag11_hs","SampleTag12_hs"),]$Cell_Index

P7_seurat$cell_hashing_cond <- NA
P7_seurat@meta.data <- P7_seurat@meta.data %>%
  mutate(cell_hashing_cond = case_when(
    rownames(P7_seurat@meta.data) %in% P7_multiplets ~ "Doublet",
    rownames(P7_seurat@meta.data) %in% P7_singlets ~ "Singlet",
    TRUE ~ NA_character_))

## extract only Doublet and Singlet to remove Undetermined cells from the objetc  
Idents(P7_seurat) <- "cell_hashing_cond"
P7_seurat <- subset(P7_seurat, idents = c("Doublet","Singlet"))

## Merge samples
obj <- merge(P1_seurat, y = c(P2_seurat, P3_seurat,P4_seurat,P5_seurat,P6_seurat,P7_seurat),
             add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7"))
obj <- JoinLayers(obj)

### Add mitochondrial percentage per cell 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

### Save object
saveRDS(obj, file = "/scratch/khandl/technical/seurat_objects/Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7.rds")

##### Load R object 
obj <- readRDS(file = "/scratch/khandl/technical/seurat_objects/Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
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
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

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
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" # endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:27)
new.cluster.ids <- c("Neutrophils","Eosinophils","T","PCs", "lowQ", "Macrophages",
                     "lowQ","Monocytes_DC","PCs", "TAMs","Epithelial",
                     "Neutrophils","Endothelial","Mast","Fibroblasts","B",
                     "lowQ","PCs","Mixed","T","lowQ",
                     "Mixed","Mast","Mixed","Mixed","Mixed","Mixed","Mixed")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")

### Subcluster Monocytes_DC
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Monocytes_DC",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Monocytes_DC_0","Monocytes_DC_1","Monocytes_DC_2","Monocytes_DC_3"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
FeaturePlot(sub_celltype, features = "CD1C")
FeaturePlot(sub_celltype, features = "VCAN")
# 0 = Monocytes, 1 = DCs, 2 = TAMs, 3 = Mixed

# Rename
current.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macrophages","Mast","Mixed","Monocytes_DC_0","Monocytes_DC_1","Monocytes_DC_2","Monocytes_DC_3", "Neutrophils","PCs","T","TAMs")
new.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages","Mast","Mixed","Monocytes","DCs","TAMs", "Mixed", "Neutrophils","PCs","T","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### Check annotation
Idents(subCl) <- "annotation"

DotPlot(subCl, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
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
                                 "COL1A2", # Fibroblasts
                                 "PECAM1" # Endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

VlnPlot(subCl, features = "nFeature_RNA")

##### Save object 
saveRDS(subCl, "/scratch/khandl/technical/seurat_objects/Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_annotated.rds")

