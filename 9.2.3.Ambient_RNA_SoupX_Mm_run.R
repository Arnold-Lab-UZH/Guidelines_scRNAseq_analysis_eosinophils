########## This code quantifies ambient RNA content using SoupX  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

# GSE282765 
# take forced cell determination and intronic + exonic reads
##### read R objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

## add batch label for each experiment so you can fun decontX on batches 
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids<- c("batch1","batch2","batch1","batch3")
obj$batch <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### run separate for each batch 
### batch 1 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/Tumor_phil_disseminated_phil_disseminated_wt_adjacent_colon_wt_Expression_Data_Unfiltered.st")

Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch1")

tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 

# tod and toc must have the same number of genes 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]

# Create the SoupChannel object
sc <- SoupChannel(tod = tod, toc = toc)

## add clustering information 
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)

# Estimate rho
sc = autoEstCont(sc) # rhos = 0.01 
# Clean the data
out_batch1 = adjustCounts(sc,roundToInt=TRUE)

batch1_Seurat <- CreateSeuratObject(out_batch1)

## assess the level of contamination 
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("Ighg1", "Ighg3","Ighd","Ighe","Ighm","Iglc1","Iglc2","Iglc3","Iglc4","Igkc")

# define the cell cluster which is allowed to express these genes, in our case PCs, these are not use to estimate the contamination 
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
# Ensure ordering matches SoupX cell order
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)

## generate a dataframe for the contamination fraction per cell type cluster 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

# adjacent_colon_wt
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "adjacent_colon_wt")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
# make data frame 
contamination_values <- c(10.77,1.00,6.58,3.01,0.96,1.83,4.49,3.34,1.43,80.6)
df1 <- data.frame(cell_types, contamination_values)

# disseminated_wt
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "disseminated_wt")
cell_types <- c("B","DCs","Eosinophils","Macrophages","Monocytes","Neutrophils","T","TAMs","PCs") # no mast cells
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
# make data frame 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast", "Monocytes","Neutrophils","T","TAMs","PCs") 
contamination_values <- c(11.49,0.79,13.50,2.10,0.00, 1.34,4.14,2.40,0.88,99)
df2 <- data.frame(cell_types, contamination_values)

### batch 2 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/Colon_healthy_wt_phil_Expression_Data_Unfiltered.st")

Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch2")

tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 

# tod and toc must have the same number of genes 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]

# Create the SoupChannel object
sc <- SoupChannel(tod = tod, toc = toc)

## add clustering information 
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)

# Estimate rho
sc = autoEstCont(sc) # rhos = 0.02 
# Clean the data
out_batch2 = adjustCounts(sc,roundToInt=TRUE)
batch2_Seurat <- CreateSeuratObject(out_batch2)

## assess the level of contamintation 
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("Ighg1", "Ighg3","Ighd","Ighe","Ighm","Iglc1","Iglc2","Iglc3","Iglc4","Igkc")

# define the cell cluster which is allowed to express these genes, in our case PCs, these are not use to estimate the contamination 
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
# Ensure ordering matches SoupX cell order
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)

## generate a dataframe for the contamination fraction per cell type cluster 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
# make data frame 
contamination_values <- c(9.81,2.56,10.15,2.85,3.28,5.31,6.10,4.87,4.53,82.49)
df3 <- data.frame(cell_types, contamination_values)

### batch 3
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/Tumor_wt_Expression_Data_Unfiltered.st")

Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch3")

tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 

# tod and toc must have the same number of genes 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]

# Create the SoupChannel object
sc <- SoupChannel(tod = tod, toc = toc)

## add clustering information 
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)

# Estimate rho
sc = autoEstCont(sc) # rhos = 0.02 
# Clean the data
out_batch3 = adjustCounts(sc,roundToInt=TRUE)
batch3_Seurat <- CreateSeuratObject(out_batch3)

### show the global degree of contamination based on igGenes 
# define the gene set 
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("Ighg1", "Ighg3","Ighd","Ighe","Ighm","Iglc1","Iglc2","Iglc3","Iglc4","Igkc")

# define the cell cluster which is allowed to express these genes, in our case PCs, these are not use to estimate the contamination 
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
# Ensure ordering matches SoupX cell order
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)

## generate a dataframe for the contamination fraction per cell type cluster 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
# make data frame 
contamination_values <- c(30.53,2.24,10.06,2.93,5.42,3.33,7.02,7.32,2.94,77.69)
df4 <- data.frame(cell_types, contamination_values)

##### plot degree of contamination acorss samples 
df <- rbind(df1,df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/SoupX_mM.csv")

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + scale_y_break(c(35, 75))   + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values = c("Eosinophils" = "#E22F27",  "Neutrophils" = "#9518ED",  
                               "B" = "#F3E972", "T" = "#5BC7D9",   "Monocytes" = "#ADD8AB", "Mast" = "#7F7F79", "PCs" = "#B4C108",
                               "Macrophages" = "#82C341", 
                               "TAMs" = "#516D38","DCs" = "#E43794","lowQ"="#282525","Mixed" = "#E59A38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/SoupX_contamination_per_celltypeMm.svg", width = 10, height = 8, plot = p)

anova <- aov(contamination_values ~ cell_types, data = df)
summary(anova)
TukeyHSD(anova)

##### merge the counts and integrate the assay to Seurat 
merged <- merge(batch1_Seurat,c(batch2_Seurat, batch3_Seurat))
merged$condition <- obj$condition
merged$percent.mt <- PercentageFeatureSet(merged, pattern = "^mt.")
merged <- JoinLayers(merged)

##### pre-process and cluster 
merged <- NormalizeData(merged,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged, features = VariableFeatures(object =merged), npcs = 20, verbose = FALSE)

### fastMNN integration 
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$condition)
merged <- IntegrateLayers(object = merged, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                          verbose = FALSE)
ElbowPlot(merged)         
merged <- FindNeighbors(merged, reduction = "integrated.mnn", dims = 1:15, graph.name = "integrated.mnn_snn")
merged <- FindClusters(merged, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2,graph.name = "integrated.mnn_snn")
merged <- RunUMAP(merged, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(merged,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
merged <- JoinLayers(merged)

##### identify the eos cluster and annotate it 
p <- FeaturePlot(merged, features = "Ccr3", reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/SoupX_Ccr3.svg", width = 10, height = 8, plot = p)

current.cluster.ids <- c(0:20)
new.cluster.ids <- c("other","other","other","other", "other", "other",
                     "other","Eosinophils","other","other","other",
                     "other","other","other","other","other",
                     "other","other", "other","other","other")
merged$eos_cluster <- plyr::mapvalues(x = merged$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, group.by = "eos_cluster", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### transfer annotation from original object 
cell_types <- (as.data.frame(table(obj$annotation)))$Var1

cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

shared <- intersect(rownames(merged@meta.data),  rownames(obj@meta.data)) #30114
additional <- rownames(merged@meta.data)[!rownames(merged@meta.data) %in% shared] #0

merged$annotation <- NA
merged@meta.data <- merged@meta.data %>%
  mutate(annotation = case_when(
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(merged$annotation)
DimPlot(merged, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### save object 
saveRDS(merged, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_SoupX.rds")
