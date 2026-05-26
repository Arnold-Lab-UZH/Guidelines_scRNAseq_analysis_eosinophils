########### This code compares the usage of  Malat1 to identify empty droplets  ##########
### Datasets used: GSE282765; BM and blood Mm healthy and CRC 

##### Set up environment 
setwd("/data/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.4.Functions_MALAT1_threshold_empty_droplets_Clarke_Bader_etal.R")

##### Run for each experiment/cartridge individually 
counts_data <- data_to_sparse_matrix_unfiltered("/scratch/khandl/Technical3/Mm_BM_blood_tumor_healthy_Expression_Data_Unfiltered.st")
counts_data_seurat <- CreateSeuratObject(counts_data)
counts_data_seurat <- NormalizeData(counts_data_seurat)

# Get normalized data matrix
norm_data <- GetAssayData(counts_data_seurat, layer = "data")
# Extract MALAT1 expression (gene names are case-sensitive!)
malat1_expr <- as.numeric(norm_data["Malat1", ])

## Calculate threshold 
threshold <- define_malat1_threshold_ggplot2(malat1_expr) #2.45795

## Apply the threshold and define real and empty cells 
df <- FetchData(counts_data_seurat, vars = c("Malat1"),layer = "data")
colnames(df) <- c("norm_expression")
df$threshold_malat1 <- ifelse(df$norm_expression > threshold, "MALAT1_real", "MALAT1_empty")
df <- df[df$threshold_malat1 %in% "MALAT1_real",]

## Identify cell barcodes from NAT and tuomr 
ST_calls <- read.csv("/scratch/khandl/Technical3/Mm_BM_blood_tumor_healthy_Sample_Tag_Calls.csv",skip = 7)
blood_healthy_calls <- (ST_calls[ST_calls$Sample_Name %in% "blood_ctrl",])$Cell_Index
blood_tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "blood_tumor",])$Cell_Index
bm_healthy_calls <- (ST_calls[ST_calls$Sample_Name %in% "BM_ctrl",])$Cell_Index
bm_tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "BM_tumor",])$Cell_Index

## Extract the right barcodes per sample from the MALAT1 analysis df 
df_blood_healthy <- df[rownames(df) %in% blood_healthy_calls,]
df_blood_tumor <- df[rownames(df) %in% blood_tumor_calls,]
df_bm_healthy <- df[rownames(df) %in% bm_healthy_calls,]
df_bm_tumor <- df[rownames(df) %in% bm_tumor_calls,]

counts_data_blood_healthy <- counts_data[,rownames(df_blood_healthy)]
counts_data_blood_tumor <- counts_data[,rownames(df_blood_tumor)]
counts_data_bm_healthy <- counts_data[,rownames(df_bm_healthy)]
counts_data_bm_tumor <- counts_data[,rownames(df_bm_tumor)]

## Generate Seurat objects
blood_wt <- create_seurat_Mm_datafrom_sparse_matrix(counts_data_blood_healthy,  project = "blood_wt", condition = "blood_wt",3,200)
blood_tumor <- create_seurat_Mm_datafrom_sparse_matrix(counts_data_blood_tumor,  project = "blood_tumor", condition = "blood_tumor",3,200)
bm_wt <- create_seurat_Mm_datafrom_sparse_matrix(counts_data_bm_healthy, project = "bm_wt", condition = "bm_wt",3,200)
bm_tumor <- create_seurat_Mm_datafrom_sparse_matrix(counts_data_bm_tumor,  project = "bm_tumor", condition = "bm_tumor",3,200)

##### Merge all Seurat objects 
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

### Add conditions to metadata 
tumor$cell_determination <- "MALAT1"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

### Apply mitochondrial cutoff 
tumor <- subset(tumor, subset = percent.mt < 25)

##### Load annotated object from BD forced pipeline 
obj_forced <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Clustering 
### Pre-processing 
obj <- tumor
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)         
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)
obj_MALAT1 <- obj

##### Transfer of annotation based on matching cell IDs 
### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?","Basophils","DCs","EoP", "Eosinophils","GMPs", "lowQ", "Macrophages","HSCs",
                         "Monocytes","Neutrophils","ProMono","ProNeutro", "T")
new.cluster.ids <- c("Undefined","Basophils","DCs","EoP", "Eosinophils","GMPs", "lowQ", "Macrophages","HSCs",
                     "Monocytes","Neutrophils","ProMono","ProNeutro", "T")
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

## Check if the cell type labels are matchin between analyses 
shared <- intersect(rownames(obj_MALAT1@meta.data),  rownames(obj_forced@meta.data)) #24133
additional_ones_in_malat1 <- rownames(obj_MALAT1@meta.data)[!rownames(obj_MALAT1@meta.data) %in% shared] #0

## Add shared and unique barcode identity to Meta data 
obj_MALAT1$annotation <- NA
obj_MALAT1@meta.data <- obj_MALAT1@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Basophils ~ "Basophils",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$EoP ~ "EoP",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$GMPs ~ "GMPs",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$HSCs ~ "HSCs",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$ProMono ~ "ProMono",
    rownames(obj_MALAT1@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_MALAT1@meta.data)  %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj_MALAT1$annotation)
DimPlot(obj_MALAT1, group.by = "annotation", label = TRUE)

##### save object 
saveRDS(obj_MALAT1, "/scratch/khandl/technical/seurat_objects/Mm_blood_bm_healthy_tumor_MALAT1_determination_with_intronic_reads_annotated.rds")



