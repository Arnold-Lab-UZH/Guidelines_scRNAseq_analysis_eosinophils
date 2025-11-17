########### This code quantifies ambient RNA content using DecontX default (cell clusters) ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

# GSE282765 
# take forced cell determination and intronic + exonic reads

##### read R objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

## add batch label for each experiment so you can fun decontX on batches 
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids<- c("batch1","batch2","batch1","batch3")
obj$batch <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### DecontX without supplying the empty droplet background, estimation of contamination based on profiels of other cell clusters in the filtered data = DecontX_default
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce <- decontX(sce,z=sce$annotation, batch = sce$batch)

##### plot contamination in feature Plot 
# add contamination result to meta data 
obj$decontX_contamination <- colData(sce)$decontX_contamination
Idents(obj) <- "condition"
sub <- subset(obj, idents = "tumor_wt")
p <- FeaturePlot(sub, features = c("decontX_contamination")) +scale_color_gradientn( colours = c('darkblue',"yellow" ,'darkred'))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_feature_plot_Mm.svg", width = 10, height = 8, plot = p)

##### plot level of contamination based on DecontX_default 
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
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/decontX_mM.csv")

p <- ggplot(df, aes(x = reorder(celltype,median_decontX, FUN = median), y =  median_decontX, fill = celltype)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values = c("Eosinophils" = "#E22F27",  "Neutrophils" = "#9518ED",  
                               "B" = "#F3E972", "T" = "#5BC7D9",   "Monocytes" = "#ADD8AB", "Mast" = "#7F7F79", "PCs" = "#B4C108",
                               "Macrophages" = "#82C341", 
                               "TAMs" = "#516D38","DCs" = "#E43794","lowQ"="#282525","Mixed" = "#E59A38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/DecontX_contamination_per_celltypeMm.svg", width = 10, height = 8, plot = p)

## statistical test --> one way ANOVA 
anova <- aov(median_decontX ~ celltype, data = df)
summary(anova)
TukeyHSD(anova)

##### convert back to a Seurat object and cluster/annotate based on original annotaiton 
## round because Seurat needs integers 
merged <- CreateSeuratObject(round(decontXcounts(sce)))
merged$condition <- obj$condition
merged$percent.mt <- PercentageFeatureSet(merged, pattern = "^mt.")

##### pre-process and cluster 
# Default Assay is RNA, but RNA assay is the decontX 

## fix error from ScaleData: 0%Error in qr.resid(qr = qr, y = data.expr[x, ]) : 'qr' and 'y' must have the same number of rows
summary(is.na(merged@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])) #2746 NA in percent.mt
keep_cells <- complete.cases(merged@meta.data[, c("percent.mt", "nCount_RNA", "nFeature_RNA")])
merged <- subset(merged, cells = Cells(merged)[keep_cells])

merged <- NormalizeData(merged,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged, features = VariableFeatures(object =merged), npcs = 20, verbose = FALSE)

### fastMNN integration 
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$condition)
merged <- IntegrateLayers(object = merged, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                          verbose = FALSE)
ElbowPlot(merged)         
merged <- FindNeighbors(merged, reduction = "integrated.mnn", dims = 1:15, graph.name = "integrated.mnn_snn")
merged <- FindClusters(merged, resolution = 0.8, cluster.name = "mnn.clusters", algorithm = 2,graph.name = "integrated.mnn_snn")
merged <- RunUMAP(merged, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(merged,reduction = "umap.mnn",group.by = "mnn.clusters",raster=TRUE, label = TRUE, label.size = 8)
merged <- JoinLayers(merged)

##### identify the eos cluster and annotate it 
p <- FeaturePlot(merged, features = "Ccr3", reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/decontX_Ccr3.svg", width = 10, height = 8, plot = p)

current.cluster.ids <- c(0:20)
new.cluster.ids <- c("other","other","other","other", "other", "other",
                     "other","Eosinophils","other","other","other",
                     "other","other","other","other","other",
                     "other","other", "other","other","other")
merged$eos_cluster <- plyr::mapvalues(x = merged$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, group.by = "eos_cluster", label = TRUE,raster=FALSE,reduction = "umap.mnn")

##### transfer annotation from original object 
cell_types <- (as.data.frame(table(obj$annotation)))$Var1

cell_ids_per_annotation_list <- list()
for (i in cell_types){
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  barcodes <- rownames(sub@meta.data)
  cell_ids_per_annotation_list[[i]] <- barcodes
}

shared <- intersect(rownames(merged@meta.data),  rownames(obj@meta.data)) #30114
additional <- rownames(merged@meta.data)[!rownames(merged@meta.data) %in% shared] #0

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

### save object 
saveRDS(merged, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_decontX.rds")
