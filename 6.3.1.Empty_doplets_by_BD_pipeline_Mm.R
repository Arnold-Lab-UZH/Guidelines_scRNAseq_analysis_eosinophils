########### This code compares forced and automatic BD pipeline ##########
### Datasets used: GSE282765; BM and blood Mm healthy and CRC 

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.2.Functions_Seurat_integration.R"))

##### Load annotated object from BD forced pipeline 
obj_forced <- readRDS(file.path(seurat_objects_dir, "Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))

##### Seurat object generation from BD automatic 
### automatic cell determination - intronic and exonic reads 
blood_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_automatic_intron_exon_dir, "Automatic_cell_determination_intronic_and_exonic_Mm_blood_healthy_ST11_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

blood_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_automatic_intron_exon_dir, "Automatic_cell_determination_intronic_and_exonic_Mm_blood_AKPS_tumor_ST12_Expression_Data.st"), 
  project = "blood_tumor", condition = "blood_tumor",3,200)

bm_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_automatic_intron_exon_dir, "Automatic_cell_determination_intronic_and_exonic_Mm_bm_healthy_ST09_Expression_Data.st"), 
  project = "bm_wt", condition = "bm_wt",3,200)

bm_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path(raw_data_GSE282765_Mm_bm_blood_automatic_intron_exon_dir, "Automatic_cell_determination_intronic_and_exonic_Mm_bm_AKPS_tumor_ST10_Expression_Data.st"), 
  project = "bm_tumor", condition = "bm_tumor",3,200)

### Merge samples
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt-")

### Add conditions to metadata 
tumor$cell_determination <- "automatic"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(tumor, file = file.path(seurat_objects_dir, "Automatic_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds"))

##### Load BD automatic object 
obj_autom <- readRDS( file.path(data_dir, "Automatic_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds"))

##### Clustering 
### Pre-processing 
obj <- obj_autom
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- subset(obj, subset = percent.mt < 25)
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
obj_autom <- obj

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
shared <- intersect(rownames(obj_autom@meta.data),  rownames(obj_forced@meta.data)) #27735
additional_ones_in_autom <- rownames(obj_autom@meta.data)[!rownames(obj_autom@meta.data) %in% shared] #0

## Add cell type labels based on cell IDs from BD forced 
obj_autom$annotation <- NA
obj_autom@meta.data <- obj_autom@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Basophils ~ "Basophils",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$EoP ~ "EoP",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$GMPs ~ "GMPs",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$HSCs ~ "HSCs",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$ProMono ~ "ProMono",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj_autom$annotation)
DimPlot(obj_autom, group.by = "annotation", label = TRUE)

##### save object 
saveRDS(obj_autom, file.path(seurat_objects_dir, "Mm_blood_bm_healthy_tumor_automatic_cell_determination_with_intronic_reads_annotated.rds"))
