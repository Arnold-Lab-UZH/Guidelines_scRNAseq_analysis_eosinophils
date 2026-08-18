######### This code identifies and analyses doublets based on DoubletFinder  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load annotated data 
obj <- readRDS(file.path(seurat_objects_dir,"Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated.rds"))

##### Run for each experiment separately 
### Exp1 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp1")

### Identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.18

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.18_2717
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.18_2717
sub_exp1 <- sub

### Exp2 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp2")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.28

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.28_2440
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.28_2440
sub_exp2 <- sub

### Exp3
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp3")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.18

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.18_2131
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.18_2131
sub_exp3 <- sub

### Exp4
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp4")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.001

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.001, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.001_1162
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.001_1162
sub_exp4<- sub

### Exp5
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp5")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.29

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.29_2242
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.29_2242
sub_exp5 <- sub

### Exp7 = patient 6
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp7")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.24

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.24_3156
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.24_3156
sub_exp7 <- sub

### Exp8 = patient7
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp8")

### identify pk 
sweep.res.list <- paramSweep(sub, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
print(optimal_pK) #0.15

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(sub$annotation)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(homotypic.prop*nrow(sub@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
sub <- doubletFinder(sub, PCs = 1:10, pN = 0.25, pK = 0.15, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
sub$DoubletFinder_score <- sub$pANN_0.25_0.15_3474
sub$DoubletFinder_class <- sub$DF.classifications_0.25_0.15_3474
sub_exp8 <- sub

### merge all 
obj <- merge(sub_exp1,y= c(sub_exp2,sub_exp3,sub_exp4,sub_exp5,sub_exp7,sub_exp8))
obj <- JoinLayers(obj)

##### pre-process and cluster 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### fastMNN integration 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)         
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "annotation", label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

Idents(obj) <- "annotation"
sub <- subset(obj, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","Macrophages","Mast","Monocytes","Neutrophils","PCs","T","TAMs"))
Idents(sub) <- "experiment"
sub <- subset(sub, idents = "Exp2")
p <- FeaturePlot(sub, features = "DoubletFinder_score", reduction = "umap.mnn",cols = c("white","darkred") )
ggsave(file.path(doublet_plots_dir, "DoubletFinder_umap.svg"), width = 8, height = 8, plot = p)

##### calculate percentage of doublet rate 
# experiment 7 = h6/P6, experiment 8 = H7/P7
experimet_ids <- c("Exp1","Exp2","Exp3","Exp4", "Exp5","Exp7", "Exp8")
for (i in experimet_ids) {
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = c(i))
  all_cells <- length(rownames(sub@meta.data))
  Idents(sub) <- "DoubletFinder_class"
  sub2 <- subset(sub, idents = "Doublet")
  doublet_cells <- length(rownames(sub2@meta.data))
  multiplet_rate <- (doublet_cells *100)/all_cells
  print(multiplet_rate)
}

sample <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
multiplet_rate_per_sample <- c(12.90491,12.34505,13.30877,12.18924,9.668794,11.0134,11.04189)
df <- data.frame(sample, multiplet_rate_per_sample)
df$method <- "DoubletFinder"
write.csv(df,file.path(doublet_tables_dir, "DoubletFinder_doublet_rate_corr.csv"))

##### extract doublets and deconvolute 
Idents(obj) <- "DoubletFinder_class"
sub <- subset(obj, idents = "Doublet")

### deconvolute multiplets with SCDC
## load and format data, should be a dataframe with genes as rownames and colums are the samples/cells 
# extract raw counts 
df <- as.matrix(sub@assays$RNA@layers$counts)
rownames(df) <- rownames(sub@assays$RNA)
colnames(df) <- colnames(sub)
df <- as.data.frame(df)

### generate reference (need NAT, tumor)
hs <- readRDS(file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))

Idents(hs) <- "annotation"
hs <- subset(hs, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts", "Macrophages",
                            "Mast","Monocytes","Neutrophils","PCs","T","TAMs"))

reference <- hs
#select 200 cells per annotated cluster 
Idents(reference) <- reference$annotation
reference <- subset(reference, cells = WhichCells(reference, downsample = 200))

## Select genes for deconvolution 
reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(reference, only.pos = TRUE, logfc.threshold = 0.25,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")

# Filter for genes that are also present in the bulk RNAseq data 
markers_sc <- markers_sc[markers_sc$gene %in% rownames(df), ]

# Select top 20 genes per cluster, select top by first p-value, then absolute diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))

### deconvolution
scRNAseq_matrix <- as.matrix(reference[["RNA"]]$counts)
# extract m_feats 
scRNAseq_matrix <-scRNAseq_matrix[rownames(scRNAseq_matrix) %in% m_feats,]

# extract m_feats from bulk data 
df_matrix <- as.matrix(df)
df_matrix <- df_matrix[rownames(df_matrix) %in% m_feats,]

# remove all columns/cells where less than 10 genes are experessed 
keep_cols <- colSums(df_matrix > 0) >= 10
df_matrix_filtered <- df_matrix[, keep_cols] #14813 lost 

# run deconvolution 
eset_SC <- ExpressionSet(assayData = scRNAseq_matrix, 
                         phenoData = AnnotatedDataFrame(reference@meta.data))
eset_ST <- ExpressionSet(df_matrix_filtered)

# remova a problematic cell and run again 
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "annotation",
                                     ct.sub = as.character(unique(eset_SC$annotation)))

deconvolution_crc_df <- as.data.frame(deconvolution_crc$prop.est.mvw)
write.csv(deconvolution_crc_df,file.path(doublet_tables_dir, "DoubletFinder_wo_doublets_deconvolution_result_corr.csv"))

##### save Seurat object 
saveRDS(obj, file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated_DoubletFinder.rds"))
