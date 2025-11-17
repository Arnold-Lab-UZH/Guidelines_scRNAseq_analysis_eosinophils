########## This code quantifies ambient RNA content  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

# GSE282765 and GSE182001
# take forced cell determination and intronic + exonic reads

##### read R objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

current.cluster.ids <- c("?","B","DCs", "Eosinophils", "lowQ", "Macrophages","Mast","Mixed",
                         "Monocytes","Neutrophils","PCs","TAMs", "T")
new.cluster.ids <- c("Undefined","B","DCs", "Eosinophils", "lowQ", "Macrophages","Mast","Mixed",
                     "Monocytes","Neutrophils","PCs","TAMs", "T")
obj$annotation <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

## add batch label for each experiment so you can fun decontX on batches 
current.cluster.ids <- c("adjacent_colon_wt","adult_colon_wt","disseminated_wt","tumor_wt")
new.cluster.ids<- c("batch1","batch2","batch1","batch3")
obj$batch <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### run for each batch 
### batch 1 
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch1")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype for both samples within batch 
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
# make data frame 
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
# make data frame 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast", "Monocytes","Neutrophils","T","TAMs","PCs") 
contamination_values <- c(0.0196,0.02282,0.01806,0.01435,0.000,0.01262,0.00554,0.0035,0.00745,0.33971)
df2 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch1_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 2 
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch2")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
# make data frame 
contamination_values <- c(0.01384,0.0339,0.02238,0.01425,0.00619,0.01761,0.01377,0.0123,0.01585,0.27712)
df3 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch2_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 3
Idents(obj) <- "batch"
sub <- subset(obj, idents = "batch3")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","T","TAMs","PCs")

for(i in cell_types) {
  Idents(sub) <- "annotation"
  sub2 <- subset(sub, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
# make data frame 
contamination_values <- c(0.0167,0.03253,0.02031,0.01976,0.00802,0.01606,0.05175,0.01388,0.01584,0.37653)
df4 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch3_corrected = ContaminationCorrection(sub,rownames(GCGs))

##### Merge samples
merged <- merge(batch1_corrected, y = c(batch2_corrected, batch3_corrected))
DefaultAssay(merged) <- "Corrected"
merged <- JoinLayers(merged)

##### pre-process and cluster 
DefaultAssay(merged) <- "Corrected"
merged <- NormalizeData(merged,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged, features = VariableFeatures(object =merged), npcs = 20, verbose = FALSE)

### fastMNN integration 
merged[["Corrected"]] <- split(merged[["Corrected"]], f = merged$condition)
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
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_Ccr3.svg", width = 10, height = 8, plot = p)

DefaultAssay(obj) <- "RNA"
p <- FeaturePlot(obj, features = "Ccr3", reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/RNA_Ccr3.svg", width = 10, height = 8, plot = p)

current.cluster.ids <- c(0:21)
new.cluster.ids <- c("other","other","other","other", "other", "other",
                     "other","other","other","other","Eosinophils",
                     "other","other","other","other","other",
                     "other","other", "other","other","other","other")
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

##### plot degree of contamination acorss samples 
df <- rbind(df1,df2)
df <- rbind(df, df3)
df <- rbind(df, df4)
write.csv(df,"/scratch/khandl/technical/figures/Ambient_RNA/scCDC_mM.csv")

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + scale_y_break(c(0.06, 0.25))   + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values = c("Eosinophils" = "#E22F27",  "Neutrophils" = "#9518ED",  
                               "B" = "#F3E972", "T" = "#5BC7D9",   "Monocytes" = "#ADD8AB", "Mast" = "#7F7F79", "PCs" = "#B4C108",
                               "Macrophages" = "#82C341", 
                               "TAMs" = "#516D38","DCs" = "#E43794","lowQ"="#282525","Mixed" = "#E59A38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_contamination_per_celltypeMm.svg", width = 10, height = 8, plot = p)

anova <- aov(contamination_values ~ cell_types, data = df)
summary(anova)
TukeyHSD(anova)

### save object 
saveRDS(merged, "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_scCDC.rds")

