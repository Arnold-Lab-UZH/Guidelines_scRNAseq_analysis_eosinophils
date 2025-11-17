######### This code compares the mesians of number of features genes of annotated cell types across datasets  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

##### human data 
### extract cell IDs from multiplets 
P1 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp1_Sample_Tag_Calls.csv",skip = 7)
P2 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp2_Sample_Tag_Calls.csv",skip = 7)
P3 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp3_Sample_Tag_Calls.csv",skip = 7)
P4 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp4_Sample_Tag_Calls.csv",skip = 7)
P5 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp5_Sample_Tag_Calls.csv",skip = 7)
P6 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp6_Sample_Tag_Calls.csv",skip = 7)
P7 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp7_Sample_Tag_Calls.csv",skip = 7)
#P8 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp_PB_healthy_Sample_Tag_Calls.csv",skip = 7)

P1 <- P1[P1$Sample_Tag %in% "Multiplet",]$Cell_Index
P2 <- P2[P2$Sample_Tag %in% "Multiplet",]$Cell_Index
P3 <- P3[P3$Sample_Tag %in% "Multiplet",]$Cell_Index
P4 <- P4[P4$Sample_Tag %in% "Multiplet",]$Cell_Index
P5 <- P5[P5$Sample_Tag %in% "Multiplet",]$Cell_Index
P6 <- P6[P6$Sample_Tag %in% "Multiplet",]$Cell_Index
P7 <- P7[P7$Sample_Tag %in% "Multiplet",]$Cell_Index
#P8 <- P8[P8$Sample_Tag %in% "Multiplet",]$Cell_Index

### from unfiltered count matrix extract the multiplets and then generat Seurat object with 200 and 3 cutoff 
## P1 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P1_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P1]
P1_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P1",3,200,  "P1_multiplet","multiplet","Exp1","patient")

## P2 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P2_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P2]
P2_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P2",3,200,  "P2_multiplet","multiplet","Exp2","patient")

## P3 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P3_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P3]
P3_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P3",3,200,  "P3_multiplet","multiplet","Exp3","patient")

## P4 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P4_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P4]
P4_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P4",3,200,  "P4_multiplet","multiplet","Exp4","patient")

## P5 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P5_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P5]
P5_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P5",3,200,  "P5_multiplet","multiplet","Exp5","patient")

## P6 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P6_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P6]
P6_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P6",3,200,  "P6_multiplet","multiplet","Exp7","patient")

## P7 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P7_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data <- counts_data[,colnames(counts_data) %in% P7]
P7_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P7",3,200,  "P7_multiplet","multiplet","Exp8","patient")

## Merge samples
obj <- merge(P1_seurat, y = c(P2_seurat, P3_seurat,P4_seurat,P5_seurat,P6_seurat,P7_seurat),
             add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7"))
obj <- JoinLayers(obj)

obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = percent.mt < 25)
length(obj@active.ident) #4943 mutliplets 

saveRDS(obj,"/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_multiplets_from_cell_hashing.rds")

## cluster and look at marker gene expression 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15, reduction.name = "umap")
DimPlot(obj,reduction = "umap", label = TRUE, label.size = 8)
FeaturePlot(obj, features = c("CLC","SYNE1","ICOS"))

### deconvolute multiplets with SCDC
## load and format data, should be a dataframe with genes as rownames and colums are the samples/cells 
# extract raw counts 
df <- as.matrix(obj@assays$RNA@layers$counts)
rownames(df) <- rownames(obj@assays$RNA)
colnames(df) <- colnames(obj)
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
write.csv(deconvolution_crc_df,"/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_deconvolution_result.csv")

### add the deconvolution result to the object 
obj <- subset(obj, cells = colnames(df_matrix_filtered))
# Make sure cell order matches
all(Cells(obj) == rownames(deconvolution_crc_df))  # should return TRUE
# Add as metadata
obj <- AddMetaData(obj, metadata = deconvolution_crc_df[, c("Endothelial","Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                            "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

##### plot the number of features across different predictions of Eos 
deconvolution_crc_df_eos_cat1 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils ==0 &deconvolution_crc_df$Eosinophils <= 0.25 ,]
deconvolution_crc_df_eos_cat2 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.25 &deconvolution_crc_df$Eosinophils <= 0.5 ,]
deconvolution_crc_df_eos_cat3 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.5 &deconvolution_crc_df$Eosinophils <= 0.75 ,]
deconvolution_crc_df_eos_cat4 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.75 ,]

## annotate the ranges 
obj$eos_scdc_prop_range <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(eos_scdc_prop_range = case_when(
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat1) ~ "cat0_25",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat2) ~ "cat25_50",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat3) ~ "cat50_75",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat4) ~ "cat75_100",
    TRUE ~ NA_character_))

## add single Eos 
all_cells <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(all_cells) <- "annotation"
singlets <- subset(all_cells, idents = c("Eosinophils"))
singlets$eos_scdc_prop_range <- singlets$annotation
obj2 <- merge(obj, singlets)

## plot 
Idents(obj2) <- "eos_scdc_prop_range"
sub <- subset(obj2, idents = c("cat0_25","cat25_50","cat50_75","cat75_100","Eosinophils"))

p <- VlnPlot(sub, features= "nFeature_RNA",pt.size = 0, cols = c("#E5D6D6","#E2ACAC" ,"#E06870","#C61825" ,"#7C0A12" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() 
print(p)
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_eos_doublet_score_nfeatures.svg", width = 8, height = 8, plot = p)

p <- DotPlot(sub, features = c("CLC","CCR3","ALOX15","ADGRE1","DACH1", "CD3E","ICOS","CD4", "C1QC","C1QB","VCAN","FN1","FCGR3B"),dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/technical/figures/Doublet/Eos_hashing_doulbets_DotPlot.svg", width = 10, height = 6, plot = p)

##### plot the number of features across different predictions of Mac 
deconvolution_crc_df_Mac_cat1 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages ==0 &deconvolution_crc_df$Macrophages <= 0.25 ,]
deconvolution_crc_df_Mac_cat2 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.25 &deconvolution_crc_df$Macrophages <= 0.5 ,]
deconvolution_crc_df_Mac_cat3 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.5 &deconvolution_crc_df$Macrophages <= 0.75 ,]
deconvolution_crc_df_Mac_cat4 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.75 ,]

## annotate the ranges 
obj$scdc_prop_rangeMac <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(scdc_prop_rangeMac = case_when(
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat1) ~ "cat0_25",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat2) ~ "cat25_50",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat3) ~ "cat50_75",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat4) ~ "cat75_100",
    TRUE ~ NA_character_))

## plot 
singlets <- subset(all_cells, idents = c("Macrophages"))
singlets$scdc_prop_rangeMac <- singlets$annotation
obj2 <- merge(obj, singlets)

## plot 
Idents(obj2) <- "scdc_prop_rangeMac"
sub <- subset(obj2, idents = c("cat0_25","cat25_50","cat50_75","cat75_100","Macrophages"))

p <- VlnPlot(sub, features= "nFeature_RNA",pt.size = 0, cols = c("#D1F2DE","#7DE2A6" ,"#1F9B51","#087505" ,"#224421" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() 
print(p)
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_Macs_doublet_score_nfeatures.svg", width = 8, height = 8, plot = p)

##### extract category 25-50 
Idents(obj) <- "eos_scdc_prop_range"
sub <- subset(obj, idents = c("cat25_50","cat50_75"))

### visualise markers 
p <- FeaturePlot(sub, features = c("CLC","ICOS","C1QC","FCGR3B"))
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_doulbet_marker.svg", width = 10, height = 8, plot = p)

### visualize deconvolution score to verify 
p <- FeaturePlot(sub, features = "Neutrophils",cols = c("lightgrey", "#5D0EAF"))
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_Neutrophils_score.svg", width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "Macrophages",cols = c("lightgrey", "#2D4F08"))
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_Mac_score.svg", width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "Eosinophils",cols = c("lightgrey", "#E23128"))
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_Eosinophils_score.svg", width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "T",cols = c("lightgrey", "#0C52A5"))
ggsave("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_T_score.svg", width = 8, height = 8, plot = p)

##### how many eos doublets are expected based on hashing 
df <- read.csv("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_deconvolution_result.csv")
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

### how many cells per sample haver larger 25 and lower 75 = doublets 
all_h1 <- length(rownames(df_h1))
doublets_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.25 & df_h1$Eosinophils <= 0.75 ,]))
homo_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.75  ,]))
other_h1 <- length(rownames(df_h1[df_h1$Eosinophils ==0 &df_h1$Eosinophils <= 0.25  ,]))
doublets_h1 <- (doublets_h1*100)/all_h1
homo_h1 <- (homo_h1*100)/all_h1
other_h1 <- (other_h1*100)/all_h1
all_h1 <- (all_h1*100)/all_h1
h1_df <- data.frame(c("all","doublets","homo","other"), c(all_h1,doublets_h1,homo_h1,other_h1))
colnames(h1_df) <- c("identity","values")
h1_df$exp <- "Exp1"

all_h2 <- length(rownames(df_h2))
doublets_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.25 & df_h2$Eosinophils <= 0.75 ,]))
homo_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.75  ,]))
other_h2 <- length(rownames(df_h2[df_h2$Eosinophils ==0 &df_h2$Eosinophils <= 0.25  ,]))
doublets_h2 <- (doublets_h2*100)/all_h2
homo_h2 <- (homo_h2*100)/all_h2
other_h2 <- (other_h2*100)/all_h2
all_h2 <- (all_h2*100)/all_h2
h2_df <- data.frame(c("all","doublets","homo","other"), c(all_h2,doublets_h2,homo_h2,other_h2))
colnames(h2_df) <- c("identity","values")
h2_df$exp <- "Exp2"

all_h3 <- length(rownames(df_h3))
doublets_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.25 & df_h3$Eosinophils <= 0.75 ,]))
homo_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.75  ,]))
other_h3 <- length(rownames(df_h3[df_h3$Eosinophils ==0 &df_h3$Eosinophils <= 0.25  ,]))
doublets_h3 <- (doublets_h3*100)/all_h3
homo_h3 <- (homo_h3*100)/all_h3
other_h3 <- (other_h3*100)/all_h3
all_h3 <- (all_h3*100)/all_h3
h3_df <- data.frame(c("all","doublets","homo","other"), c(all_h3,doublets_h3,homo_h3,other_h3))
colnames(h3_df) <- c("identity","values")
h3_df$exp <- "Exp3"

all_h4 <- length(rownames(df_h4))
doublets_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.25 & df_h4$Eosinophils <= 0.75 ,]))
homo_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.75  ,]))
other_h4 <- length(rownames(df_h1[df_h4$Eosinophils ==0 &df_h4$Eosinophils <= 0.25  ,]))
doublets_h4 <- (doublets_h4*100)/all_h4
homo_h4 <- (homo_h4*100)/all_h4
other_h4 <- (other_h4*100)/all_h4
all_h4 <- (all_h4*100)/all_h4
h4_df <- data.frame(c("all","doublets","homo","other"), c(all_h4,doublets_h4,homo_h4,other_h4))
colnames(h4_df) <- c("identity","values")
h4_df$exp <- "Exp4"

all_h5 <- length(rownames(df_h5))
doublets_h5 <- length(rownames(df_h1[df_h5$Eosinophils >0.25 & df_h5$Eosinophils <= 0.75 ,]))
homo_h5 <- length(rownames(df_h5[df_h5$Eosinophils >0.75  ,]))
other_h5 <- length(rownames(df_h5[df_h5$Eosinophils ==0 &df_h5$Eosinophils <= 0.25  ,]))
doublets_h5 <- (doublets_h5*100)/all_h5
homo_h5 <- (homo_h5*100)/all_h5
other_h5 <- (other_h5*100)/all_h5
all_h5 <- (all_h5*100)/all_h5
h5_df <- data.frame(c("all","doublets","homo","other"), c(all_h5,doublets_h5,homo_h5,other_h5))
colnames(h5_df) <- c("identity","values")
h5_df$exp <- "Exp5"

all_h6 <- length(rownames(df_h6))
doublets_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.25 & df_h6$Eosinophils <= 0.75 ,]))
homo_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.75  ,]))
other_h6 <- length(rownames(df_h6[df_h6$Eosinophils ==0 &df_h6$Eosinophils <= 0.25  ,]))
doublets_h6 <- (doublets_h6*100)/all_h6
homo_h6 <- (homo_h6*100)/all_h6
other_h6 <- (other_h6*100)/all_h6
all_h6 <- (all_h6*100)/all_h6
h6_df <- data.frame(c("all","doublets","homo","other"), c(all_h6,doublets_h6,homo_h6,other_h6))
colnames(h6_df) <- c("identity","values")
h6_df$exp <- "Exp6"

all_h7 <- length(rownames(df_h7))
doublets_h7 <- length(rownames(df_h6[df_h7$Eosinophils >0.25 & df_h7$Eosinophils <= 0.75 ,]))
homo_h7 <- length(rownames(df_h7[df_h7$Eosinophils >0.75  ,]))
other_h7 <- length(rownames(df_h6[df_h7$Eosinophils ==0 &df_h7$Eosinophils <= 0.25  ,]))
doublets_h7 <- (doublets_h7*100)/all_h7
homo_h7 <- (homo_h7*100)/all_h7
other_h7 <- (other_h7*100)/all_h7
all_h7 <- (all_h7*100)/all_h7
h7_df <- data.frame(c("all","doublets","homo","other"), c(all_h7,doublets_h7,homo_h7,other_h7))
colnames(h7_df) <- c("identity","values")
h7_df$exp <- "Exp7"

df <- rbind(h1_df,h2_df)
df <- rbind(df,h3_df)
df <- rbind(df,h4_df)
df <- rbind(df,h5_df)
df <- rbind(df,h6_df)
df <- rbind(df,h7_df)

p <- ggplot(df, aes(x = reorder(identity,values, FUN = median), y =  values, fill = identity)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,100) +
  scale_fill_manual(values =  c( "doublets" = "#DBC20B", "homo" = "#E22F27","other" = "#08991C", 
                                 "all" = "#1EC9D6"))
ggsave("/scratch/khandl/technical/figures/Doublet/Eos_doublet_rate_hashing.svg", width = 8, height = 8, plot = p)

