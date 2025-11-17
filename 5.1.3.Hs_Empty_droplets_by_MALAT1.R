########### This code compares the usage of  MALAT1 to identify empty droplets  ##########

##### Set up environment 
setwd("/data/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.4.Functions_MALAT1_threshold_empty_droplets_Clarke_Bader_etal.R")

##### plot Malat1 expression in different cell types 
### Homo sapiens 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

#### plot expression of MALAT1 in histogram and visualize the distribution 
gene <- "MALAT1"
obj <- NormalizeData(obj)
df <- FetchData(obj, vars = c(gene,"annotation"),layer = "data")
colnames(df) <- c("norm_expression","annotation")

p <- ggplot(df, aes(x = norm_expression, fill = annotation == "Eosinophils")) +
  geom_histogram(binwidth = 0.1, color = "black", alpha = 0.7) +
  theme_minimal() +
  xlab("Normalized Expression of Malat1") +
  ylab("Number of Cells") +
  ggtitle("Distribution of Malat1 Expression") +
  scale_fill_manual(values = c("grey", "#E22F27"), labels = c("Other", "cell_lower200"))
ggsave("/scratch/khandl/technical/figures/Empty_droplets_by_Malat1/Hs_Malat1_histogram_Eos_highlighted.svg", width = 12, height = 8, plot = p)

##### run for each experiment individually 
### P1 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P1_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.597286

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp1_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "contol",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P1_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

### P2 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P2_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.475979

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp2_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "h2_tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "h2_tissue_ctrl",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P2_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P2",3,200, "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

### P3 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P3_tumor_NAT_blood_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.247958

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp3_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "h3_tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "h3_control_CD45",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P3_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P3",3,200, "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

### P4
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P4_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.545198

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp4_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "control_CD45",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P4_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P4",3,200, "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

### P5
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P5_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.454798

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp5_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor_CD45",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "control_CD45",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P5_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P5",3,200, "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

### P6 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P6_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.518954

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp6_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "P6_tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "P6_control",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P6_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P6",3,200,  "P6_tumor","tumor","Exp6","patient")
P6_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P6",3,200, "P6_tissue_ctrl","tissue_ctrl","Exp6","patient")

### P7
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical_count_matrices/P7_tumor_NAT_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["MALAT1", ])

## calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.296447

## add threshold to cell barcodes 
df <- FetchData(counts_data_seurat, vars = c("MALAT1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp7_Sample_Tag_Calls.csv",skip = 7)
tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "tumor",])$Cell_Index
NAT_calls <- (ST_calls[ST_calls$Sample_Name %in% "tissue_ctrl",])$Cell_Index

## extract the right barcodes per sample from the EmptyDrops df 
df_tumor <- df[rownames(df) %in% tumor_calls,]
df_NAT <- df[rownames(df) %in% NAT_calls,]

counts_data_tumor <- counts_data[,rownames(df_tumor)]
counts_data_NAT <- counts_data[,rownames(df_NAT)]

## make Seurat objects 
P7_tumor <- create_seurat_Hs_data_from_sparse_matrix(counts_data_tumor, "P7",3,200,  "P7_tumor","tumor","Exp7","patient")
P7_control <- create_seurat_Hs_data_from_sparse_matrix(counts_data_NAT, "P7",3,200, "P7_tissue_ctrl","tissue_ctrl","Exp7","patient")

## Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

patients$cell_determination <- "MALAT1"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

patients <- subset(patients, subset = percent.mt < 25)

##### add cell type label from forced BD annotation
obj_forced <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### pre-process and cluster 
obj <- patients
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)         
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)
obj_MALAT1 <- obj

##### add annotation based on cell id match 
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
obj_forced$annotation <- plyr::mapvalues(x = obj_forced$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_forced$annotation)))$Var1

cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_forced) <- "annotation"
  sub <- subset(obj_forced, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

shared <- intersect(rownames(obj_MALAT1@meta.data),  rownames(obj_forced@meta.data)) #121590
additional_ones_in_malat1 <- rownames(obj_MALAT1@meta.data)[!rownames(obj_MALAT1@meta.data) %in% shared] #1

obj_MALAT1$annotation <- NA
obj_MALAT1@meta.data <- obj_MALAT1@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(obj_MALAT1$annotation)
DimPlot(obj_MALAT1, group.by = "annotation", label = TRUE)

## merge objects 
# remove NA from obj_MALAT1
Idents(obj_MALAT1) <- "annotation"
obj_MALAT1 <- subset(obj_MALAT1, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                                                    "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

### save object 
saveRDS(obj_MALAT1, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_MALAT1_determination_with_intronic_reads_annotated.rds")



