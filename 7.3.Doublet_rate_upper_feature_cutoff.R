######### This code compares the mesians of number of features genes of annotated cell types across datasets  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.6.Functions_Doublet_detection.R")

##### Hs NAT/tumor data 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### define a treshold as the upper 99th percentile of nFeature_RNA in FeatureScatter plto 
# show the example of experiment 2 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp2")
nf <- sub$nFeature_RNA
upper_q99 <- quantile(nf, 0.99)
p <- FeatureScatter(sub, "nCount_RNA", "nFeature_RNA", group.by="orig.ident", pt.size=.5)+   geom_point(color = "#0C8EEF") + 
  geom_hline(yintercept = upper_q99, color="darkred", linetype="dashed")
ggsave("/scratch/khandl/technical/figures/Doublet/FeatureScatter_example.svg", width = 10, height = 8, plot = p)

##### annotate the doublets based on the individual cutoff 
exp1 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp1") 
exp2 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp2") 
exp3 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp3") 
exp4 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp4") 
exp5 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp5") 
exp6 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp7") 
exp7 <- annotate_doublets_based_on_upper99th_percentile(obj, "Exp8") 

### combine in the original object 
Idents(exp1) <- "doublet_singlet"
exp1_doublet <- subset(exp1, idents = "doublet")
exp1_doublet <- rownames(exp1_doublet@meta.data)
exp1_singlet <- subset(exp1, idents = "singlet")
exp1_singlet <- rownames(exp1_singlet@meta.data)

Idents(exp2) <- "doublet_singlet"
exp2_doublet <- subset(exp2, idents = "doublet")
exp2_doublet <- rownames(exp2_doublet@meta.data)
exp2_singlet <- subset(exp2, idents = "singlet")
exp2_singlet <- rownames(exp2_singlet@meta.data)

Idents(exp3) <- "doublet_singlet"
exp3_doublet <- subset(exp3, idents = "doublet")
exp3_doublet <- rownames(exp3_doublet@meta.data)
exp3_singlet <- subset(exp3, idents = "singlet")
exp3_singlet <- rownames(exp3_singlet@meta.data)

Idents(exp4) <- "doublet_singlet"
exp4_doublet <- subset(exp4, idents = "doublet")
exp4_doublet <- rownames(exp4_doublet@meta.data)
exp4_singlet <- subset(exp4, idents = "singlet")
exp4_singlet <- rownames(exp4_singlet@meta.data)

Idents(exp5) <- "doublet_singlet"
exp5_doublet <- subset(exp5, idents = "doublet")
exp5_doublet <- rownames(exp5_doublet@meta.data)
exp5_singlet <- subset(exp5, idents = "singlet")
exp5_singlet <- rownames(exp5_singlet@meta.data)

Idents(exp6) <- "doublet_singlet"
exp6_doublet <- subset(exp6, idents = "doublet")
exp6_doublet <- rownames(exp6_doublet@meta.data)
exp6_singlet <- subset(exp6, idents = "singlet")
exp6_singlet <- rownames(exp6_singlet@meta.data)

Idents(exp7) <- "doublet_singlet"
exp7_doublet <- subset(exp7, idents = "doublet")
exp7_doublet <- rownames(exp7_doublet@meta.data)
exp7_singlet <- subset(exp7, idents = "singlet")
exp7_singlet <- rownames(exp7_singlet@meta.data)

obj$doublet_singlet <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(doublet_singlet = case_when(
    rownames(obj@meta.data) %in% c(exp1_doublet,exp2_doublet,exp3_doublet,exp4_doublet,exp5_doublet,exp6_doublet,exp7_doublet) ~ "doublet",
    rownames(obj@meta.data) %in% c(exp1_singlet,exp2_singlet,exp3_singlet,exp4_singlet,exp5_singlet,exp6_singlet,exp7_singlet) ~ "singlet",
    TRUE ~ NA_character_))

##### calculate doublet rate per experiment 
experimet_ids <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
for (i in experimet_ids) {
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = c(i))
  all_cells <- length(rownames(sub@meta.data))
  Idents(sub) <- "doublet_singlet"
  sub2 <- subset(sub, idents = "doublet")
  doublet_cells <- length(rownames(sub2@meta.data))
  multiplet_rate <- (doublet_cells *100)/all_cells
  print(multiplet_rate)
}

sample <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
multiplet_rate_per_sample <- c(1.002381,1.000702,1.005818,1.001431,1.002297,1.001263,1.000266)
df <- data.frame(sample, multiplet_rate_per_sample)
df$method <- "nFeature_upper_cutoff"
write.csv(df,"/scratch/khandl/technical/figures/Doublet/nFeature_upper_cutoff_doublet_rate.csv")

##### extract doublets and deconvolute 
Idents(obj) <- "doublet_singlet"
sub <- subset(obj, idents = "doublet")

### deconvolute multiplets with SCDC
## load and format data, should be a dataframe with genes as rownames and colums are the samples/cells 
# extract raw counts 
df <- as.matrix(sub@assays$RNA@layers$counts)
rownames(df) <- rownames(sub@assays$RNA)
colnames(df) <- colnames(sub)
df <- as.data.frame(df)

### generate reference (need NAT, tumor)
hs <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

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
write.csv(deconvolution_crc_df,"/scratch/khandl/technical/figures/Doublet/upperFeature_cutoff_doublets_deconvolution_result.csv")

##### save Seurat object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated_gene_count_cutoff.rds")

