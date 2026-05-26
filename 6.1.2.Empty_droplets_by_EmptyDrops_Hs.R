########### This code compares forced and the EmptyDrops pipeline ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor;

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### run emptyDrops by cartridge/experiment -> 0.01 FDR cutoff 
### P1 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P1_tumor_NAT_blood_Expression_Data_Unfiltered.st")
set.seed(100)
# Test ambient means that cells with lower than 100 UMI counts are also tested against the ambient profile 
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]

## Assign cell barcodes based on Sample Tag Calls
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P1_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "contol",])$Cell_Index

## Extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## Generate Seurat objects 
P1_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

### P2
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P2_tumor_NAT_blood_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P2_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "h2_tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "h2_tissue_ctrl",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P2_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P2",3,200, "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

### P3
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P3_tumor_NAT_blood_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P3_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "h3_tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "h3_control_CD45",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P3_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P3",3,200, "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

### P4
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P4_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P4_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "control_CD45",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P4_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P4",3,200, "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

### P5
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P5_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P5_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "control_CD45",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P5_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P5",3,200, "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

### P6
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P6_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P6_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "P6_tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "P6_control",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P6_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P6",3,200, "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

### P7
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P7_tumor_NAT_Expression_Data_Unfiltered.st")
set.seed(100)
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_P7_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "tissue_ctrl",])$Cell_Index
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]
counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]
P7_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P7",3,200, "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

###### Plot an example scatter plot highlighting eosinophils as empty drops and real 
## Rerun  
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)

### Extract Eos to plot them in the EmptyDrops plot 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

# Extract eosinophil specific cell barcodes 
Idents(obj) <- "annotation"
eos <- subset(obj, idents = "Eosinophils")

Idents(eos) <- "experiment"
sub <- subset(eos, idents = "Exp8")
eos_bc <- rownames(sub@meta.data)
eos_bc <- substr(eos_bc, 5, nchar(eos_bc))
df_TRUE <- df[df$is.cell == TRUE,]
df_FALSE <- df[df$is.cell == FALSE,]
matching_bc_real <- intersect(rownames(df_TRUE), eos_bc)
matching_bc_empty <- intersect(rownames(df_FALSE), eos_bc)

# Logical vectors for each group
is_real <- rownames(df) %in% matching_bc_real
is_empty <- rownames(df) %in% matching_bc_empty

# Define colors for other cells
colors <- ifelse(df$is.cell, "blue", "black")

# Assign colors for the two special groups
colors[is_real] <- "red"
colors[is_empty] <- "yellow"

# Plot
plot(df$Total, -df$LogProb, col = colors, 
     xlab = "Total UMI count", ylab = "-Log Probability")
points(df$Total[is_real], -df$LogProb[is_real], col = "red", cex = 1.2, pch = 16)
points(df$Total[is_empty], -df$LogProb[is_empty], col = "yellow", cex = 1.2,pch = 16)

##### Merge all Seurat objects 
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

### Add mitochondrial percentage per cell 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

### Add conditions to metadata 
patients$cell_determination <- "EmptyDrops"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

### Apply mitochondrial cutoff 
patients <- subset(patients, subset = percent.mt < 25)

##### Load annotated object from BD forced pipeline 
obj_forced <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Clustering 
### Pre-processing 
obj <- patients
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
obj_emptyDrops <- obj

##### Transfer of annotation based on matching cell IDs 
### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
obj_forced$annotation <- plyr::mapvalues(x = obj_forced$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_forced$annotation)))$Var1

## Loop through cell type labels and extract cell IDs for each label 
cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_forced) <- "annotation"
  sub <- subset(obj_forced, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

## Check if the cell type labels are matching between analyses 
shared <- intersect(rownames(obj_emptyDrops@meta.data),  rownames(obj_forced@meta.data)) #122704
additional_ones_in_emptyDrops <- rownames(obj_emptyDrops@meta.data)[!rownames(obj_emptyDrops@meta.data) %in% shared] #1

## Add shared and unique barcode identity to Meta data 
obj_emptyDrops$annotation <- NA
obj_emptyDrops@meta.data <- obj_emptyDrops@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(obj_emptyDrops$annotation)
DimPlot(obj_emptyDrops, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_forced, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj_emptyDrops, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

# Remove NA 
Idents(obj_emptyDrops) <- "annotation"
obj_emptyDrops <- subset(obj_emptyDrops, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                                                    "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

##### save object 
saveRDS(obj_emptyDrops, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_emptyDrops_determination_with_intronic_reads_annotated.rds")



