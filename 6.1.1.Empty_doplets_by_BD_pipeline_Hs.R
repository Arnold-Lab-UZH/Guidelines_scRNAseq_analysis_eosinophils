########### This code compares forced and automatic BD pipeline ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor;

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.2.Functions_Seurat_integration.R"))

##### Load annotated object from BD forced pipeline 
obj_forced <- readRDS(file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))

##### Seurat object generation from BD automatic 
### Automatic cell determination - intronic and exonic reads 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P1_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                    "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P1_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P2_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P2_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P3_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P3_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P4_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P4_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P5_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P5_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P6_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P6_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir, "P7_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path(raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir, "P7_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

### Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

### Add mitochondrial percentage per cell 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

### Add conditions to metadata 
patients$cell_determination <- "automatic"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

### Save object
saveRDS(patients, file = file.path(data_dir, "Automatic_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds"))

##### Load BD automatic object 
obj_autom <- readRDS( file.path(data_dir, "Automatic_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds"))

### Apply mitochondrial cutoff 
obj_autom <- subset(obj_autom, subset = percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- obj_autom
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$experiment)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
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
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 8)
obj <- JoinLayers(obj)
obj_autom <- obj

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
shared <- intersect(rownames(obj_autom@meta.data),  rownames(obj_forced@meta.data)) #78669
additional_ones_in_autom <- rownames(obj_autom@meta.data)[!rownames(obj_autom@meta.data) %in% shared] #673

## Add shared and unique barcode identity to Meta data 
obj_autom$shared_unique <- NA
obj_autom@meta.data <- obj_autom@meta.data %>%
  mutate(shared_unique = case_when(
    rownames(obj_autom@meta.data) %in% shared ~ "shared_BC",
    rownames(obj_autom@meta.data)  %in% additional_ones_in_autom ~ "unique_BD_autom",
    TRUE ~ NA_character_))
table(obj_autom$shared_unique)
DimPlot(obj_autom, group.by = "shared_unique")

## Add cell type labels based on cell IDs from BD forced 
obj_autom$annotation <- NA
obj_autom@meta.data <- obj_autom@meta.data %>%
  mutate(annotation = case_when(
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Endothelial ~ "Endothelial",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Epithelial ~ "Epithelial",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Fibroblasts ~ "Fibroblasts",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(obj_autom@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(obj_autom@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(obj_autom$annotation)
DimPlot(obj_autom, group.by = "annotation", label = TRUE)

## Subset additional barcodes from BD automatic, cluster and annotate 
Idents(obj_autom) <- "shared_unique"
sub <- subset(obj_autom, idents = "unique_BD_autom")
DimPlot(sub, group.by = "mnn.clusters", label =TRUE)

Idents(sub) <- "mnn.clusters"
DotPlot(sub, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1","FFAR2", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", #N eutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" #Endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 
VlnPlot(sub, features = "nFeature_RNA")

# Renam clusters 
current.cluster.ids <- c(0:28)
new.cluster.ids <- c("PCs","T","Neutrophils", "Macrophages","lowQ", "TAMs", 
                     "Monocytes","Epithelial","Fibroblasts","PCs", "Fibroblasts",
                     "Mixed","Mast","T","Mixed","Mixed",
                     "Mixed","Neutrophils","Mixed","Mixed","Mixed",
                     "Mixed","Mixed","Mixed","Mast","Mixed",
                     "Mixed","Mixed","Mixed")
sub$annotation <- plyr::mapvalues(x = sub$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
p1 <- DimPlot(sub, group.by = "annotation", label = TRUE,reduction = "umap.mnn")
p2 <- DimPlot(obj_autom, group.by = "annotation", label = TRUE,reduction = "umap.mnn")
p1 + p2

## Merge objects from the shared cell barcodes and the unique barcodes from BD automatic 
# Remove NA from obj_autom
Idents(obj_autom) <- "annotation"
obj_autom <- subset(obj_autom, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ","Macrophages","Mast","Mixed",
                                          "Monocytes","Neutrophils","PCs","T","TAMs","Undefined"))

obj_autom <- merge(obj_autom, sub)
obj_autom <- JoinLayers(obj_autom)

##### save object 
saveRDS(obj_autom, file.path(seurat_objects_dir, "Hs_tumor_NAT_automatic_cell_determination_with_intronic_reads_annotated.rds"))

