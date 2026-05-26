########### This code compares forced BD pipeline with two different gene mapping strategies (exons + introns and exons only )  ##########
### Datasets used: GSE282765; Mm colon healthy, tumor, disseminated NAT 

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Load annotated object from BD forced pipeline with intron and exon mapping 
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Seurat object generation from BD forced exons only  
### Forced cell determination exonic reads only 
tumor_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_tumor_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_disseminated_ST07_Expression_Data.st"), 
  project = "disseminated_wt", condition = "disseminated_wt",3,200)

adult_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_colon_ST04_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

adjacent_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_NAT_ST03_Expression_Data.st"), 
  project = "adjacent_colon_wt", condition = "adjacent_colon_wt",3,200)

### Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt,adjacent_colon_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt","adjacent_colon_wt"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell  
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

### Add conditions to metadata 
tumor$cell_determination <- "forced"
tumor$reads <- "exonic_only"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

### Save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Mm_CRC_AKPS_tumor_NAT_diss.rds")

##### Load exons only object 
obj <- readRDS("/scratch/khandl/4.Technical/Forced_cell_determination_exonic_reads_only_Mm_CRC_AKPS_tumor_NAT_diss.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
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
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label = TRUE, label.size = 8)
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
obj <- JoinLayers(obj)

##### Transfer of annotation based on matching cell IDs 
### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?", "B","DCs","Eosinophils",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs","Eosinophils",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
obj_reference$annotation <- plyr::mapvalues(x = obj_reference$annotation, from = current.cluster.ids, to = new.cluster.ids)

cell_types <- (as.data.frame(table(obj_reference$annotation)))$Var1

## Loop through cell type labels and extract cell IDs for each label 
cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj_reference) <- "annotation"
  sub <- subset(obj_reference, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

## Check if the cell type labels are matching between analyses 
shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #49269
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #0

## Add cell type labels based on cell IDs from BD forced 
obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(obj$annotation)
DimPlot(obj, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_reference, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

# Remove NA
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B","DCs","Eosinophils","lowQ","Macrophages","Mast","Mixed",
                              "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

##### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_exons_only_annotation.rds")



