######### This code identifies and analyses doublets based on sample tag cell hashing results ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.2.Functions_Seurat_integration.R"))

##### Extract multiplets from Unfiltered matrices based on Sample Tag Calls 
### Load Sample Tag Calls and extract Multiplets Cell IDs 
P1 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P1_Sample_Tag_Calls.csv"),skip = 7)
P2 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P2_Sample_Tag_Calls.csv"),skip = 7)
P3 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P3_Sample_Tag_Calls.csv"),skip = 7)
P4 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P4_Sample_Tag_Calls.csv"),skip = 7)
P5 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P5_Sample_Tag_Calls.csv"),skip = 7)
P6 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P6_Sample_Tag_Calls.csv"),skip = 7)
P7 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P7_Sample_Tag_Calls.csv"),skip = 7)

P1 <- P1[P1$Sample_Tag %in% "Multiplet",]$Cell_Index
P2 <- P2[P2$Sample_Tag %in% "Multiplet",]$Cell_Index
P3 <- P3[P3$Sample_Tag %in% "Multiplet",]$Cell_Index
P4 <- P4[P4$Sample_Tag %in% "Multiplet",]$Cell_Index
P5 <- P5[P5$Sample_Tag %in% "Multiplet",]$Cell_Index
P6 <- P6[P6$Sample_Tag %in% "Multiplet",]$Cell_Index
P7 <- P7[P7$Sample_Tag %in% "Multiplet",]$Cell_Index

### Load unfiltered matrices and extract Multiplets Cell IDs and generate Seurat objects 
## P1 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P1_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P1]
P1_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P1",3,200,  "P1_multiplet","multiplet","Exp1","patient")

## P2 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P2_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P2]
P2_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P2",3,200,  "P2_multiplet","multiplet","Exp2","patient")

## P3 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P3_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P3]
P3_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P3",3,200,  "P3_multiplet","multiplet","Exp3","patient")

## P4 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P4_tumor_NAT_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P4]
P4_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P4",3,200,  "P4_multiplet","multiplet","Exp4","patient")

## P5 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P5_tumor_NAT_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P5]
P5_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P5",3,200,  "P5_multiplet","multiplet","Exp5","patient")

## P6 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P6_tumor_NAT_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P6]
P6_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P6",3,200,  "P6_multiplet","multiplet","Exp7","patient")

## P7 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"P7_tumor_NAT_Expression_Data_Unfiltered.st"))
counts_data <- counts_data[,colnames(counts_data) %in% P7]
P7_seurat <- create_seurat_Hs_data_from_sparse_matrix(counts_data, "P7",3,200,  "P7_multiplet","multiplet","Exp8","patient")

## Merge samples
obj <- merge(P1_seurat, y = c(P2_seurat, P3_seurat,P4_seurat,P5_seurat,P6_seurat,P7_seurat),
             add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7"))
obj <- JoinLayers(obj)

## Apply mitochondrial percentage 
obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

## Apply mitochondiral cutoff 
obj <- subset(obj, subset = percent.mt < 25)
length(obj@active.ident) # 4943 Multiplets 

## Save Multiplets Seurat object 
saveRDS(obj,file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_multiplets_from_cell_hashing.rds"))

##### Clustering 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,  margin = 1,assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15, reduction.name = "umap")
DimPlot(obj,reduction = "umap", label = TRUE, label.size = 8)
FeaturePlot(obj, features = c("CLC","SYNE1","ICOS"))

##### Deconvolution using SCDC 
# Extract raw counts from object with doublets 
df <- as.matrix(obj@assays$RNA@layers$counts)
rownames(df) <- rownames(obj@assays$RNA)
colnames(df) <- colnames(obj)
df <- as.data.frame(df)

### Generate reference (need NAT, tumor)
reference <- readRDS(file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))

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
# extract m_feats 
scRNAseq_matrix <-scRNAseq_matrix[rownames(scRNAseq_matrix) %in% m_feats,]

# Extract m_feats from bulk data 
df_matrix <- as.matrix(df)
df_matrix <- df_matrix[rownames(df_matrix) %in% m_feats,]

# Remove all columns/cells where less than 10 genes are expressed 
keep_cols <- colSums(df_matrix > 0) >= 10
df_matrix_filtered <- df_matrix[, keep_cols]

# Run deconvolution 
eset_SC <- ExpressionSet(assayData = scRNAseq_matrix, 
                         phenoData = AnnotatedDataFrame(reference@meta.data))
eset_ST <- ExpressionSet(df_matrix_filtered)
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "annotation",
                                     ct.sub = as.character(unique(eset_SC$annotation)))

## Save the results for each doublet in a datframe 
deconvolution_crc_df <- as.data.frame(deconvolution_crc$prop.est.mvw)
write.csv(deconvolution_crc_df,file.path(doublet_tables_dir, "cell_hashing_doublets_deconvolution_result.csv"))

### Add the deconvolution result to the doublet containing Seurat object  
obj <- subset(obj, cells = colnames(df_matrix_filtered))
# Make sure cell order matches
all(Cells(obj) == rownames(deconvolution_crc_df))  # should return TRUE
# Add as metadata
obj <- AddMetaData(obj, metadata = deconvolution_crc_df[, c("Endothelial","Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                            "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

##### Plot the number of genes across different Eos estimates 
deconvolution_crc_df_eos_cat1 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils ==0 &deconvolution_crc_df$Eosinophils <= 0.25 ,]
deconvolution_crc_df_eos_cat2 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.25 &deconvolution_crc_df$Eosinophils <= 0.5 ,]
deconvolution_crc_df_eos_cat3 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.5 &deconvolution_crc_df$Eosinophils <= 0.75 ,]
deconvolution_crc_df_eos_cat4 <- deconvolution_crc_df[deconvolution_crc_df$Eosinophils >0.75 ,]

## Annotate the ranges 
obj$eos_scdc_prop_range <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(eos_scdc_prop_range = case_when(
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat1) ~ "cat0_25",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat2) ~ "cat25_50",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat3) ~ "cat50_75",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_eos_cat4) ~ "cat75_100",
    TRUE ~ NA_character_))

## Add single Eos as a control  
all_cells <- readRDS(file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))

Idents(all_cells) <- "annotation"
singlets <- subset(all_cells, idents = c("Eosinophils"))
singlets$eos_scdc_prop_range <- singlets$annotation
obj2 <- merge(obj, singlets)

## Plot counts and marker gene expression 
Idents(obj2) <- "eos_scdc_prop_range"
sub <- subset(obj2, idents = c("cat0_25","cat25_50","cat50_75","cat75_100","Eosinophils"))

p <- VlnPlot(sub, features= "nFeature_RNA",pt.size = 0, cols = c("#E5D6D6","#E2ACAC" ,"#E06870","#C61825" ,"#7C0A12" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() 
print(p)
ggsave(file.path(doublet_plots_dir, "cell_hashing_eos_doublet_score_nfeatures.svg"), width = 8, height = 8, plot = p)

p <- DotPlot(sub, features = c("CLC","CCR3","ALOX15","ADGRE1","DACH1", "CD3E","ICOS","CD4", "C1QC","C1QB","VCAN","FN1","FCGR3B"),dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave(file.path(doublet_plots_dir, "Eos_hashing_doulbets_DotPlot.svg"), width = 10, height = 6, plot = p)

##### Plot the number of features across different predictions of Macrophages 
deconvolution_crc_df_Mac_cat1 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages ==0 &deconvolution_crc_df$Macrophages <= 0.25 ,]
deconvolution_crc_df_Mac_cat2 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.25 &deconvolution_crc_df$Macrophages <= 0.5 ,]
deconvolution_crc_df_Mac_cat3 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.5 &deconvolution_crc_df$Macrophages <= 0.75 ,]
deconvolution_crc_df_Mac_cat4 <- deconvolution_crc_df[deconvolution_crc_df$Macrophages >0.75 ,]

## Annotate the ranges 
obj$scdc_prop_rangeMac <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(scdc_prop_rangeMac = case_when(
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat1) ~ "cat0_25",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat2) ~ "cat25_50",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat3) ~ "cat50_75",
    rownames(obj@meta.data) %in% rownames(deconvolution_crc_df_Mac_cat4) ~ "cat75_100",
    TRUE ~ NA_character_))

singlets <- subset(all_cells, idents = c("Macrophages"))
singlets$scdc_prop_rangeMac <- singlets$annotation
obj2 <- merge(obj, singlets)

## Plot 
Idents(obj2) <- "scdc_prop_rangeMac"
sub <- subset(obj2, idents = c("cat0_25","cat25_50","cat50_75","cat75_100","Macrophages"))

p <- VlnPlot(sub, features= "nFeature_RNA",pt.size = 0, cols = c("#D1F2DE","#7DE2A6" ,"#1F9B51","#087505" ,"#224421" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() 
print(p)
ggsave(file.path(doublet_plots_dir, "cell_hashing_Macs_doublet_score_nfeatures.svg"), width = 8, height = 8, plot = p)

##### Extract category 25-50 
Idents(obj) <- "eos_scdc_prop_range"
sub <- subset(obj, idents = c("cat25_50","cat50_75"))

### visualise markers 
p <- FeaturePlot(sub, features = c("CLC","ICOS","C1QC","FCGR3B"))
ggsave(file.path(doublet_plots_dir, "cell_hashing_doulbet_marker.svg"), width = 10, height = 8, plot = p)

### Visualize deconvolution score to verify 
p <- FeaturePlot(sub, features = "Neutrophils",cols = c("lightgrey", "#5D0EAF"))
ggsave(file.path(doublet_plots_dir, "cell_hashing_doublets_Neutrophils_score.svg"), width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "Macrophages",cols = c("lightgrey", "#2D4F08"))
ggsave(file.path(doublet_plots_dir, "cell_hashing_doublets_Mac_score.svg"), width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "Eosinophils",cols = c("lightgrey", "#E23128"))
ggsave(file.path(doublet_plots_dir, "cell_hashing_doublets_Eosinophils_score.svg"), width = 8, height = 8, plot = p)

p <- FeaturePlot(sub, features = "T",cols = c("lightgrey", "#0C52A5"))
ggsave(file.path(doublet_plots_dir, "cell_hashing_doublets_T_score.svg"), width = 8, height = 8, plot = p)
