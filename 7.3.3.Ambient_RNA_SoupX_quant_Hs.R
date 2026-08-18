########## This code uses SoupX for ambient RNA quantification  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.2.Functions_Seurat_integration.R"))

##### Load R object 
obj <- readRDS( file.path(seurat_objects_dir,"Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))

### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                         "Monocytes","Neutrophils","PCs","TAMs", "T")
new.cluster.ids <- c("Undefined","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                     "Monocytes","Neutrophils","PCs","TAMs", "T")
obj$annotation <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

##### Run for each experiment  
### Exp1 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P1_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
  
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp1")

tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 

# Tod and toc must have the same number of genes 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]

# Create the SoupChannel object
sc <- SoupChannel(tod = tod, toc = toc)

## Add clustering information 
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)

# Estimate rho
sc = autoEstCont(sc) # rhos = 0.01 
# Clean the data
out_batch1 = adjustCounts(sc,roundToInt=TRUE)

batch1_Seurat <- CreateSeuratObject(out_batch1)

## Assess the level of contamination 
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

# Define the cell cluster which is allowed to express these genes, in our case PCs, these are not use to estimate the contamination 
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
# Ensure ordering matches SoupX cell order
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)

## Generate a dataframe for the contamination fraction per cell type cluster 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P1_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(10.64,3.46,4.33,13.32,4.3,3.05,4.16,4.26,4.53,11.09,66.03,3.13,6.56)
df1 <- data.frame(cell_types, contamination_values)

# Tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P1_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(13.84,8.82,13.01,11.38,6.00,7.50,7.91,10.76,15.14,14.28,73.01,4.6,13.95)
df2 <- data.frame(cell_types, contamination_values)

### Exp2
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P2_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp2")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
# Clean the data
out_batch2 = adjustCounts(sc,roundToInt=TRUE)
batch2_Seurat <- CreateSeuratObject(out_batch2)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P2_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(21.79,15.36,4.83,18.88,11.37,6.42,9.37,12.06,16.77,14.45,91.12,7.63,8.84)
df3 <- data.frame(cell_types, contamination_values)
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P2_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(20.59,3.17,6.17,17.93,3.19,3.08,3.45,5.89,4.47,10.81,83.22,3.71,6.54)
df4 <- data.frame(cell_types, contamination_values)

### Exp 3
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P3_tumor_NAT_blood_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp3")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
out_batch3 = adjustCounts(sc,roundToInt=TRUE)
batch3_Seurat <- CreateSeuratObject(out_batch3)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P3_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(34.97,11.06,10.36,26.11,10.07,11.51,9.05,23.40,10.01,54.88,99,17.33,20.38)
df5 <- data.frame(cell_types, contamination_values)
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P3_tumor")
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","TAMs", "T")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(67.71,6.68,6.87,27.77,10.81,10.63,8.79,14.25,11.97,24.87,99,6.57,18.34)
df6 <- data.frame(cell_types, contamination_values)

### Exp 4
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P4_tumor_NAT_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp4")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
out_batch4 = adjustCounts(sc,roundToInt=TRUE)
batch4_Seurat <- CreateSeuratObject(out_batch4)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P4_tissue_ctrl")
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","TAMs", "T")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(3.08,2.22,2.54,13.45,1.92,2.74,3.11,3.95,2.49,5.26,99, 1.86,2.71)
df7 <- data.frame(cell_types, contamination_values)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P4_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(5.33,2.08,4.03,12.83,1.70,2.85,2.35,2.45, 2.89,7.23,77.71,2.01,7.45)
df8 <- data.frame(cell_types, contamination_values)

### Exp5 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P5_tumor_NAT_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp5")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
out_batch5 = adjustCounts(sc,roundToInt=TRUE)
batch5_Seurat <- CreateSeuratObject(out_batch5)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P5_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(116.67,4.10,11.86,16.35,3.20,4.47,3.32,6.64,5.89,12.36,95.97,2.74,10.33)
df9 <- data.frame(cell_types, contamination_values)
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P5_tumor")
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","TAMs", "T")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(22.93,3.42,4.44,16.06,2.42,3.17,3.31,6.47,4.49,11.46,99, 3.09,7.04)
df10 <- data.frame(cell_types, contamination_values)

### Exp7 = patient 6
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P6_tumor_NAT_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp7")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
out_batch6 = adjustCounts(sc,roundToInt=TRUE)
batch6_Seurat <- CreateSeuratObject(out_batch6)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P6_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(8.43,3.33,2.87,10.83,2.23,2.36,2.76,3.30,4.66,7.62,70.55,6.78,5.87)
df11 <- data.frame(cell_types, contamination_values)
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P6_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(10.32,2.18,2.58,9.25,1.69,1.66,2.29,2.69,2.67,6.33,75.78,2.03,4.15)
df12 <- data.frame(cell_types, contamination_values)

### Exp 8 = patient 7
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_dir, "1.GSE282765/Unfiltered/P7_tumor_NAT_Expression_Data_Unfiltered.st"))
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp8")
tod <- counts_data       # ambient profile 
toc <- sub@assays$RNA$counts  # filtered counts 
common_genes <- intersect(rownames(tod), rownames(toc))
length(common_genes)
tod <- tod[common_genes, , drop = FALSE]
toc <- toc[common_genes, , drop = FALSE]
sc <- SoupChannel(tod = tod, toc = toc)
meta <- sub@meta.data
sc <- setClusters(sc, setNames(meta$annotation, rownames(meta)))
umap_coords <- sub@reductions$umap.mnn@cell.embeddings
sc <- setDR(sc, umap_coords)
sc = autoEstCont(sc) # rhos = 0.01 
out_batch7 = adjustCounts(sc,roundToInt=TRUE)
batch7_Seurat <- CreateSeuratObject(out_batch7)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20) # Igkc is the top gene 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
clusters <- sc$metaData$cluster
names(clusters) <- rownames(sc$metaData)   # assign names = cell IDs
clusters <- clusters[colnames(sc$toc)]
useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes), cluster = FALSE)
plotMarkerMap(sc, geneSet = igGenes, useToEst = useToEst)
sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P7_tissue_ctrl")
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(22.29,12.24,12.97,21.15,0,8.40,9.13,11.85,8.97,16.01,86.49,7.37,14.90)
df13 <- data.frame(cell_types, contamination_values)
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P7_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  sub2_cell_ids <- rownames(sub2@meta.data)
  
  subset_mat <- useToEst[rownames(useToEst) %in% sub2_cell_ids, , drop = FALSE]
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = subset_mat,forceAccept = TRUE)
  print(i)
}
contamination_values <- c(15.28,8.54,10.95,20.32,10.82,9.60,8.98,12.32,8.73,15.37,76.72,7.27,11.69)
df14 <- data.frame(cell_types, contamination_values)

##### Combine all data frames and save for later 
df_list <- list(df1, df2, df3, df4, df5, df6, df7,df8, df9, df10, df11, df12, df13, df14)
df <- do.call(rbind, df_list) 
write.csv(df,file.path(ambient_rna_tables_dir, "SoupX_Hs.csv"))
