########## This code uses DecontX for ambient RNA detection and decontamination  ##########
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

## Add batch label so you can run decontX  per batch/cartridge 
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids<- c("batch1","batch2","batch1","batch3")
obj$batch <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### Run DecontX 
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce <- decontX(sce,z=sce$annotation, batch = sce$batch)

##### plot contamination in feature Plot 
### Add the contamination result to the meta.data from the Seurat object 
obj$decontX_contamination <- colData(sce)$decontX_contamination
Idents(obj) <- "condition"
sub <- subset(obj, idents = "tumor_wt")
p <- FeaturePlot(sub, features = c("decontX_contamination")) +scale_color_gradientn( colours = c('darkblue',"yellow" ,'darkred'))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_feature_plot_Mm.svg", width = 10, height = 8, plot = p)

### Generate a dataframe of decontX_contamination for each sample and cell types 
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_decontX = median(decontX_contamination, na.rm = TRUE)) %>%
    as.data.frame()
  df$condition <- i 
  df$species <- "Mm"
  df_list[[i]] <- df
}
df <- bind_rows(df_list)

df <- df[df$celltype %in% c("Eosinophils","Neutrophils","B","T","Monocytes",
                            "Mast","PCs","Macrophages","TAMs","DCs"),]
### save data frame 
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/decontX_mM.csv")

##### Convert the decontaminated count matrix back to a Seurat object 
### Round because Seurat needs integers 
merged <- CreateSeuratObject(round(decontXcounts(sce)))
merged$condition <- obj$condition

### Add mitochondrial ratio 
merged$percent.mt <- PercentageFeatureSet(merged, pattern = "^mt.")

##### Clustering 
### Pre-processing 
## Fix error from ScaleData: 0%Error in qr.resid(qr = qr, y = data.expr[x, ]) : 'qr' and 'y' must have the same number of rows
summary(is.na(merged@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])) #2746 NA in percent.mt
keep_cells <- complete.cases(merged@meta.data[, c("percent.mt", "nCount_RNA", "nFeature_RNA")])
merged <- subset(merged, cells = Cells(merged)[keep_cells])

merged[["RNA"]] <- split(merged[["RNA"]], f = merged$condition)
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
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/decontX_Ccr3.svg", width = 10, height = 8, plot = p)

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
saveRDS(merged, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_decontX.rds")
