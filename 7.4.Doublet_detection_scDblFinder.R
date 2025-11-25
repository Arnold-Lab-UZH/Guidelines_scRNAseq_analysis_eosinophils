######### This code identifies and analyses doublets based on scDblFinder  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load annotated data 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Run scDblFinder  
set.seed(100)
sce <- as.SingleCellExperiment(obj)
colData(sce)
sce <- scDblFinder(sce, clusters="annotation", sample = "experiment", dbr = 0.025)
table(truth=sce$annotation, call=sce$scDblFinder.class)

## Add the scDblFinder.score and scDblFinder.class to Seurat object 
df_dbl <- as.data.frame(colData(sce))
df_dbl$cell_barcode <- colnames(sce)
df_dbl <- df_dbl[, c("cell_barcode", "annotation", "scDblFinder.class", "scDblFinder.score")]
# Make sure cell order matches
all(Cells(obj) == df_dbl$cell_barcode)  # should return TRUE
# Add as metadata
obj <- AddMetaData(obj, metadata = df_dbl[, c("scDblFinder.class", "scDblFinder.score")])

### Plot scDblFinder score 
Idents(obj) <- "annotation"
sub <- subset(obj, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","Macrophages","Mast","Monocytes","Neutrophils","PCs","T","TAMs"))

p <- FeaturePlot(obj, features = "scDblFinder.score", reduction = "umap.mnn", cols = brewer.pal(n = 9, name = "YlOrRd") )  # gradient from low -> mid -> high )
ggsave("/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_umap0025.svg", width = 8, height = 8, plot = p)

##### Calculate percentage of doublet rate 
experimet_ids <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
for (i in experimet_ids) {
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = c(i))
  all_cells <- length(rownames(sub@meta.data))
  Idents(sub) <- "scDblFinder.class"
  sub2 <- subset(sub, idents = "doublet")
  doublet_cells <- length(rownames(sub2@meta.data))
  multiplet_rate <- (doublet_cells *100)/all_cells
  print(multiplet_rate)
}

sample <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
multiplet_rate_per_sample <- c(6.509209,6.554307,5.985603,5.612413,6.022134,6.210536,6.347697)
df <- data.frame(sample, multiplet_rate_per_sample)
df$method <- "scDblFinder"
write.csv(df,"/scratch/khandl/technical/figures/Doublet/scDblFinder_doublet_rate0025.csv")

##### Extract doublets and deconvolute 
Idents(obj) <- "scDblFinder.class"
sub <- subset(obj, idents = "doublet")

## Load and format data, should be a dataframe with genes as rownames and colums are the samples/cells 
# Extract raw counts 
df <- as.matrix(sub@assays$RNA@layers$counts)
rownames(df) <- rownames(sub@assays$RNA)
colnames(df) <- colnames(sub)
df <- as.data.frame(df)

### Generate reference (need NAT, tumor)
reference <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

Idents(reference) <- "annotation"
reference <- subset(reference, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts", "Macrophages",
                            "Mast","Monocytes","Neutrophils","PCs","T","TAMs"))

# Select 200 cells per annotated cluster 
Idents(reference) <- reference$annotation
reference <- subset(reference, cells = WhichCells(reference, downsample = 200))

## Select genes for deconvolution 
reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(reference, only.pos = TRUE, logfc.threshold = 0.25,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")

# Filter for genes that are also present in the bulk RNAseq data 
markers_sc <- markers_sc[markers_sc$gene %in% rownames(df), ]

# Select top 20 genes per cluster
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))

### Deconvolution
scRNAseq_matrix <- as.matrix(reference[["RNA"]]$counts)
scRNAseq_matrix <-scRNAseq_matrix[rownames(scRNAseq_matrix) %in% m_feats,]

# Extract m_feats from bulk data 
df_matrix <- as.matrix(df)
df_matrix <- df_matrix[rownames(df_matrix) %in% m_feats,]

# rRemove all columns/cells where less than 10 genes are expressed 
keep_cols <- colSums(df_matrix > 0) >= 10
df_matrix_filtered <- df_matrix[, keep_cols] 

eset_SC <- ExpressionSet(assayData = scRNAseq_matrix, 
                         phenoData = AnnotatedDataFrame(reference@meta.data))
eset_ST <- ExpressionSet(df_matrix_filtered)
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "annotation",
                                     ct.sub = as.character(unique(eset_SC$annotation)))

deconvolution_crc_df <- as.data.frame(deconvolution_crc$prop.est.mvw)

### Save deconvolution results 
write.csv(deconvolution_crc_df,"/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_doublets_deconvolution_result0025.csv")

##### Save Seurat object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated_scDblFinder0025.rds")
