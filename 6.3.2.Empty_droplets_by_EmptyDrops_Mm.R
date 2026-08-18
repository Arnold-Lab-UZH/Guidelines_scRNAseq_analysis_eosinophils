########### This code compares forced and the EmptyDrops pipeline ##########
### Datasets used: GSE282765; BM and blood Mm healthy and CRC 

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.2.Functions_Seurat_integration.R"))

##### run emptyDrops by cartridge/experiment -> 0.01 FDR cutoff 
### adult colon 
counts_data <- data_to_sparse_matrix_unfiltered(file.path(raw_data_GSE282765_unfiltered_dir,"Mm_BM_blood_tumor_healthy_Expression_Data_Unfiltered.st"))

set.seed(100)
# Test ambient means that cells with lower than 100 UMI counts are also tested against the ambient profile 
df <- emptyDrops(counts_data,test.ambient = TRUE)
df$is.cell <- df$FDR <= 0.01
table(df$is.cell) 
df <- as.data.frame(df)
df <- df[df$is.cell == TRUE,]

## Assign cell barcodes based on Sample Tag Calls
ST_calls <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir, "Mm_BM_blood_tumor_healthy_Sample_Tag_Calls.csv"),skip = 7)
blood_healthy_calls <- (ST_calls[ST_calls$Sample_Name %in% "blood_ctrl",])$Cell_Index
blood_tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "blood_tumor",])$Cell_Index
bm_healthy_calls <- (ST_calls[ST_calls$Sample_Name %in% "BM_ctrl",])$Cell_Index
bm_tumor_calls <- (ST_calls[ST_calls$Sample_Name %in% "BM_tumor",])$Cell_Index

## Extract the right barcodes per sample from the EmptyDrops df 
df_blood_healthy <- df[rownames(df) %in% blood_healthy_calls,]
df_blood_tumor <- df[rownames(df) %in% blood_tumor_calls,]
df_bm_healthy <- df[rownames(df) %in% bm_healthy_calls,]
df_bm_tumor <- df[rownames(df) %in% bm_tumor_calls,]

counts_data_blood_healthy <- counts_data[,rownames(df_blood_healthy)]
counts_data_blood_tumor <- counts_data[,rownames(df_blood_tumor)]
counts_data_bm_healthy <- counts_data[,rownames(df_bm_healthy)]
counts_data_bm_tumor <- counts_data[,rownames(df_bm_tumor)]

##### Generate Seurat objects 
blood_wt <- create_seurat_Mm_data_from_sparse_matrix(counts_data_blood_healthy,  project = "blood_wt", condition = "blood_wt",3,200)
blood_tumor <- create_seurat_Mm_data_from_sparse_matrix(counts_data_blood_tumor,  project = "blood_tumor", condition = "blood_tumor",3,200)
bm_wt <- create_seurat_Mm_data_from_sparse_matrix(counts_data_bm_healthy, project = "bm_wt", condition = "bm_wt",3,200)
bm_tumor <- create_seurat_Mm_data_from_sparse_matrix(counts_data_bm_tumor,  project = "bm_tumor", condition = "bm_tumor",3,200)

##### Merge all Seurat objects 
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt-")

### Add conditions to metadata 
tumor$cell_determination <- "EmptyDrops"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

### Apply mitochondrial cutoff 
tumor <- subset(tumor, subset = percent.mt < 25)

##### Load annotated object from BD forced pipeline 
obj_forced <- readRDS(file.path(seurat_objects_dir, "Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))

##### Clustering 
### Pre-processing 
obj <- tumor
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
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
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 8)
obj <- JoinLayers(obj)
obj_emptyDrops <- obj

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

## Check if the cell type labels are matching between analyses 
shared <- intersect(rownames(obj_emptyDrops@meta.data),  rownames(obj_forced@meta.data)) #30114
additional_ones_in_emptyDrops <- rownames(obj_emptyDrops@meta.data)[!rownames(obj_emptyDrops@meta.data) %in% shared] #0

## Add shared and unique barcode identity to Meta data 
obj_emptyDrops$annotation <- NA
obj_emptyDrops@meta.data <- obj_emptyDrops@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Basophils ~ "Basophils",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$EoP ~ "EoP",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$GMPs ~ "GMPs",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$HSCs ~ "HSCs",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$ProMono ~ "ProMono",
    rownames(obj_emptyDrops@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_emptyDrops@meta.data)  %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj_emptyDrops$annotation)
DimPlot(obj_emptyDrops, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_forced, group.by = "annotation", label = TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj_emptyDrops, group.by = "annotation", label = TRUE,reduction = "umap.mnn")
p1 + p2

##### save object 
saveRDS(obj_emptyDrops, file.path(seurat_objects_dir, "Mm_blood_bm_healthy_tumor_emptyDrops_determination_with_intronic_reads_annotated.rds"))



