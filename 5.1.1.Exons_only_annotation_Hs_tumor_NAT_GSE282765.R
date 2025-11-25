########### This code compares forced BD pipeline with two different gene mapping strategies (exons + introns and exons only )  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor;

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Load annotated object from BD forced pipeline with intron and exon mapping 
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Seurat object generation from BD forced exons only  
### Forced cell determination exonic reads only 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P1_Hs_tumor_forced_cell_determination_exonic_only_ST02_Expression_Data.st"), 
                                  "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P1_Hs_NAT_forced_cell_determination_exonic_only_ST01_Expression_Data.st"), 
                                    "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P2_Hs_tumor_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                  "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P2_Hs_NAT_forced_cell_determination_exonic_only_ST08_Expression_Data.st"), 
                                    "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P3_Hs_tumor_forced_cell_determination_exonic_only_ST09_Expression_Data.st"), 
                                  "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P3_Hs_NAT_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                                    "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P4_Hs_tumor_forced_cell_determination_exonic_only_ST04_Expression_Data.st"), 
                                  "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P4_Hs_NAT_forced_cell_determination_exonic_only_ST03_Expression_Data.st"), 
                                    "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P5_Hs_tumor_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                  "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P5_Hs_NAT_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                    "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P6_Hs_tumor_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                  "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P6_Hs_NAT_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                    "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P7_Hs_tumor_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                                  "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P7_Hs_NAT_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                                    "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

### Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

### Add mitochondrial percentage per cell  
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

### Add conditions to metadata 
patients$cell_determination <- "forced"
patients$reads <- "exonic_only"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

### Save object
saveRDS(patients, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Hs_NAT_tumor.rds")

##### Load exons only object 
obj <- readRDS("/scratch/khandl/4.Technical/Forced_cell_determination_exonic_reads_only_Hs_NAT_tumor.rds")

### Apply mitochondrial cutoff 
obj <- subset(obj, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)         
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
obj <- JoinLayers(obj)

##### Transfer of annotation based on matching cell IDs 
### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
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
shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #103411
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #2

## Add cell type labels based on cell IDs from BD forced 
obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
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
obj <- subset(obj, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                              "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

##### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_exons_only_annotation.rds")



