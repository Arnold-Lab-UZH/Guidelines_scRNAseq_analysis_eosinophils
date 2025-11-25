########### This code compares forced BD pipeline with two different gene mapping strategies (exons + introns and exons only )  ##########
### Datasets used: GSE282765; Hs PB healthy and CRC 

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### Load annotated object from BD forced pipeline with intron and exon mapping 
obj_reference <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")

##### Seurat object generation from BD forced exons only  
### Forced cell determination exonic reads only 
## Healthy individuals 
H1 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H1_Hs_PB_healthy_forced_cell_determination_exonic_only_ST01_Expression_Data.st"), 
                            "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H2_Hs_PB_healthy_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                            "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H3_Hs_PB_healthy_forced_cell_determination_exonic_only_ST09_Expression_Data.st"), 
                            "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H4_Hs_PB_healthy_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                            "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H5_Hs_PB_healthy_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                            "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H6_Hs_PB_healthy_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                            "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")
H7 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H7_Hs_PB_healthy_forced_cell_determination_exonic_only_ST04_Expression_Data.st"), 
                            "H7",3,200,  "H7_blood","blood_healthy","Exp1","healthy")
H8 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H8_Hs_PB_healthy_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                            "H8",3,200,  "H8_blood","blood_healthy","Exp2","healthy")
H9 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H9_Hs_PB_healthy_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                            "H9",3,200,  "H9_blood","blood_healthy","Exp3","healthy")

## CRC patient blood 
P1_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P1_Hs_PB_forced_cell_determination_exonic_only_ST03_Expression_Data.st"), 
                                  "P1",3,200,  "P1_blood","blood_patient","Exp1","patient")
P2_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P2_Hs_PB_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                  "P2",3,200,  "P2_blood","blood_patient","Exp2","patient")
P3_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P3_Hs_PB_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                                  "P3",3,200,  "P3_blood","blood_patient","Exp3","patient")
P4_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P4_Hs_PB_forced_cell_determination_exonic_only_ST02_Expression_Data.st"), 
                                  "P4",3,200,  "P4_blood","blood_patient","Exp4","patient")
P5_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P5_Hs_PB_forced_cell_determination_exonic_only_ST08_Expression_Data.st"), 
                                  "P5",3,200,  "P5_blood","blood_patient","Exp5","patient")
P6_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P6_Hs_PB_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                                  "P6",3,200,  "P6_blood","blood_patient","Exp7","patient")
P7_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P7_Hs_PB_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                                  "P7",3,200,  "P7_blood","blood_patient","Exp8","patient")

### Merge samples
blood <- merge(H1, y = c(H2,H3,H4,H5,H6,H7,H8,H9,P1_blood,P2_blood,P3_blood,P4_blood,P5_blood,P6_blood,P7_blood),
               add.cell.ids = c("b1", "b2","b3","b4", "b5","b6","b7","b8","b9","pb1", "pb2","pb3","pb4", "pb5","pb6","pb7"))
blood <- JoinLayers(blood)

### Add mitochondrial percentage per cell  
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

### Add conditions to metadata 
blood$cell_determination <- "forced"
blood$reads <- "exonic_only"
blood$species <- "Hs"
blood$technology <- "BD_Rhapsody"
blood$cell_enrichment  <- "Eosinophils"

### Save object
saveRDS(blood, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Hs_PB_healthy_and_CRC.rds")

##### Load exons only object 
obj <- readRDS("/scratch/khandl/4.Technical/Forced_cell_determination_exonic_reads_only_Hs_PB_healthy_and_CRC.rds")

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
current.cluster.ids <- c("?", "Basophils","DCs","Eosinophils","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","ProNeutro", "T")
new.cluster.ids <- c("Undefined", "Basophils","DCs","Eosinophils","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","ProNeutro", "T")
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
shared <- intersect(rownames(obj@meta.data),  rownames(obj_reference@meta.data)) #42240
additional <- rownames(obj@meta.data)[!rownames(obj@meta.data) %in% shared] #1

## Add cell type labels based on cell IDs from BD forced 
obj$annotation <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Basophils ~ "Basophils",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj@meta.data)  %in% cell_ids_per_annotation_list$ProNeutro ~ "ProNeutro",
    TRUE ~ NA_character_))
table(obj$annotation)
DimPlot(obj, group.by = "annotation", label = TRUE)

p1 <- DimPlot(obj_reference, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj, group.by = "annotation", label = TRUE,raster=TRUE,reduction = "umap.mnn")
p1 + p2

# Remove NA
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Basophils","DCs","Eosinophils","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                              "Monocytes","Neutrophils","PCs","T","ProNeutro","Undefined"))

##### save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_exons_only_annotation.rds")



