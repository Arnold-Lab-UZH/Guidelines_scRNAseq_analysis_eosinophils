########## This code uses scCDC for ambient RNA detection and decontamination  ##########
### Datasets used: GSE282765; Mm colon healthy, CRC tumor, NAT, disseminated 

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load R object 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

### Change the cell type label ? to Undefined 
current.cluster.ids <- c("?","B","DCs", "Eosinophils", "lowQ", "Macrophages","Mast","Mixed",
                         "Monocytes","Neutrophils","PCs","TAMs", "T")
new.cluster.ids <- c("Undefined","B","DCs", "Eosinophils", "lowQ", "Macrophages","Mast","Mixed",
                     "Monocytes","Neutrophils","PCs","TAMs", "T")
obj$annotation <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

## Add batch label so you can run scCDC per batch/cartridge 
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids<- c("batch1","batch2","batch1","batch3")
obj$batch <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### Run for each batch 
### Batch 1 
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch1")

## Define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## Quantify contamination per celltype for each sample within batch 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

# adjacent_colon_wt
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "adjacent_colon_wt")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
# Make data frame 
contamination_values <- c(0.0181,0.02297,0.00804,0.02399,0.00147,0.01542,0.00868,0.00552,0.01247,0.28631)
df1 <- data.frame(cell_types, contamination_values)

# disseminated_wt
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "disseminated_wt")
cell_types <- c("B","DCs","Eosinophils","Macrophages","Monocytes","Neutrophils","T","TAMs","PCs") # no mast cells
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
# Make data frame 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast", "Monocytes","Neutrophils","T","TAMs","PCs") 
contamination_values <- c(0.0196,0.02282,0.01806,0.01435,0.000,0.01262,0.00554,0.0035,0.00745,0.33971)
df2 <- data.frame(cell_types, contamination_values)

## Correct the contamination
Idents(sub) <- "annotation"
batch1_corrected = ContaminationCorrection(sub,rownames(GCGs))

### Batch 2 
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch2")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.01384,0.0339,0.02238,0.01425,0.00619,0.01761,0.01377,0.0123,0.01585,0.27712)
df3 <- data.frame(cell_types, contamination_values)

Idents(sub) <- "annotation"
batch2_corrected = ContaminationCorrection(sub,rownames(GCGs))

### Batch 3
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch3")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")
for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.0167,0.03253,0.02031,0.01976,0.00802,0.01606,0.05175,0.01388,0.01584,0.37653)
df4 <- data.frame(cell_types, contamination_values)

Idents(sub) <- "annotation"
batch3_corrected = ContaminationCorrection(sub,rownames(GCGs))

##### Merge dataframes of contamination and save dataframe for later 
df <- rbind(df1,df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/scCDC_mM.csv")

##### Merge Seurat objects 
merged <- merge(batch1_corrected, y = c(batch2_corrected, batch3_corrected))
DefaultAssay(merged) <- "Corrected"
merged <- JoinLayers(merged)

##### Clustering 
### Pre-processing 
DefaultAssay(merged) <- "Corrected"
merged[["Corrected"]] <- split(merged[["Corrected"]], f = merged$condition)
merged <- NormalizeData(merged,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged, features = VariableFeatures(object =merged), npcs = 20, verbose = FALSE)

### FastMNN integration 
merged <- IntegrateLayers(object = merged, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                          verbose = FALSE)
ElbowPlot(merged)         
merged <- FindNeighbors(merged, reduction = "integrated.mnn", dims = 1:15, graph.name = "integrated.mnn_snn")
merged <- FindClusters(merged, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2,graph.name = "integrated.mnn_snn")
merged <- RunUMAP(merged, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(merged,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
merged <- JoinLayers(merged)

##### Plot Ccr3 to identify eos in UMAP 
p <- FeaturePlot(merged, features = "Ccr3", reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_Ccr3.svg", width = 10, height = 8, plot = p)

## compare with RNA assay 
DefaultAssay(obj) <- "RNA"
p <- FeaturePlot(obj, features = "Ccr3", reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/RNA_Ccr3.svg", width = 10, height = 8, plot = p)

##### Transfer of annotation based on matching cell IDs 
cell_types <- (as.data.frame(table(obj$annotation)))$Var1

## Loop through cell type labels and extract cell IDs for each label 
cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

## Check if the cell type labels are matchin between analyses 
shared <- intersect(rownames(merged@meta.data),  rownames(obj@meta.data)) #30114
additional <- rownames(merged@meta.data)[!rownames(merged@meta.data) %in% shared] #0

## Add cell type labels based on cell IDs from BD forced 
merged$annotation <- NA
merged@meta.data <- merged@meta.data %>%
  mutate(annotation = case_when(
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Undefined ~ "Undefined",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$B ~ "B",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$DCs ~ "DCs",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Eosinophils ~ "Eosinophils",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$lowQ ~ "lowQ",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Macrophages ~ "Macrophages",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Mast ~ "Mast",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Mixed ~ "Mixed",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$Monocytes ~ "Monocytes",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$Neutrophils ~ "Neutrophils",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$PCs ~ "PCs",
    rownames(merged@meta.data) %in% cell_ids_per_annotation_list$T ~ "T",
    rownames(merged@meta.data)  %in% cell_ids_per_annotation_list$TAMs ~ "TAMs",
    TRUE ~ NA_character_))
table(merged$annotation)
DimPlot(merged, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Save object 
saveRDS(merged, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_scCDC.rds")

