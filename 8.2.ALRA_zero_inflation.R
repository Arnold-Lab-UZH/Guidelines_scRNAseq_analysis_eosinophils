########## This code compares dropout events in different cell types  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Guidelines_scRNAseq_analysis_eosinophils/Technical/1.1.Packages.R")

##### human blood dataset
obj <- readRDS( "/scratch/khandl/4.Technical/Mm_il5tg_steady_state_automatic_forced_cell_determination_with_without_intronic_reads_annotated.rds")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")

# run ALRA, creates alra assay of imputed values
obj <- RunALRA(obj)

# visualize original and imputed values
obj <- NormalizeData(obj, assay = 'RNA')
features.plot <- c( 'Siglecf')
DefaultAssay(obj) <- 'RNA'
plot1 <- FeaturePlot(obj, features.plot, ncol = 2)
DefaultAssay(obj) <- 'alra'
plot2 <- FeaturePlot(obj, features.plot, ncol = 2, cols = c('lightgrey','darkred'))
p <- CombinePlots(list(plot1, plot2), ncol = 1)
ggsave("/scratch/khandl/technical/figures/Dropouts/Eos_marker.svg", width = 10, height = 10, plot = p)

Idents(obj) <- "annotation"
eos <- WhichCells(obj, idents = c("Eosinophils"))
eop <- WhichCells(obj, idents = c("EoP"))
p <- DimPlot(obj, label=T, group.by="condition", cells.highlight= list(eos,eop), cols.highlight = c( "pink","red"), cols= "#A39F9F")
ggsave("/scratch/khandl/technical/figures/Dropouts/annotated.svg", width = 10, height = 8, plot = p)

### Compare clustering with alra and RNA assay 
Idents(obj) <- "annotation"
sub <- subset(obj, idents = c("Eosinophils","EoP"))

## RNA assay 
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
sub <- RunPCA(object = sub, features = VariableFeatures(object =sub), npcs = 20, verbose = FALSE)
ElbowPlot(sub)
sub <- FindNeighbors(sub, reduction = "pca", dims = 1:15)
sub <- FindClusters(sub, resolution = 0.2, algorithm = 2)
sub <- RunUMAP(sub, reduction = "pca", dims = 1:15, reduction.name = "umap")
DimPlot(sub,reduction = "umap",raster=FALSE, label = TRUE)
sub$RNA_clusters <- Idents(sub)
p <- DimPlot(sub, reduction = "umap", label = TRUE)
ggsave("/scratch/khandl/technical/figures/Dropouts/DimPlot_RNA_clusters.svg", width = 10, height = 8, plot = p)

p <- FeaturePlot(sub,reduction = "umap",  features = c("Epx","Ear1","Cd80","Cd274"), ncol = 2, cols = c('lightgrey','darkred'))
ggsave("/scratch/khandl/technical/figures/Dropouts/Eos_marker_RNA.svg", width = 15, height = 20, plot = p)

## annotation
markers <-  c("Mki67","Tuba1b","Epx","Prg2","Ear1","Ear2","Alox15","Aldh2","S100a9","S100a6","S100a10","Il5","Retnla","Ccl9","Il1rl1","Cd24a","Mmp9",
                          "Icosl","Il4","Tgfb1","Pirb","Rara","Cd80","Cd274","Ptgs2","Il1rn","Tnf","Siglecf")
p <- DotPlot(sub, features = markers,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/Dropouts/RNA_clus_dotPlot_markers.svg", width = 10, height = 8, plot = p)

current.cluster.ids <- c(0:8)
new.cluster.ids <- c("basal","active","?","circulating","precursor", "immature","?","?","basal")
sub$eos_subtypes <- plyr::mapvalues(x = sub$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(sub,reduction = "umap",raster=FALSE, label = TRUE, group.by = "eos_subtypes")
p <- DimPlot(sub, group.by = "eos_subtypes", reduction = "umap", label = TRUE, cols = c("#10A069","#E81818" ,"#878283","#E9EF11", "#26DFED","#E88A1A")) + ggtitle("RNA assay")
ggsave("/scratch/khandl/technical/figures/Dropouts/DimPlot_RNA_eos_subtypes.svg", width = 10, height = 8, plot = p)

## Alra assay
DefaultAssay(sub) <- "alra"
#sub <- NormalizeData(sub,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1)
sub <- FindVariableFeatures(sub)
sub <- ScaleData(sub,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
sub <- RunPCA(object = sub, features = VariableFeatures(object =sub), npcs = 20, verbose = FALSE)
ElbowPlot(sub)
sub <- FindNeighbors(sub, reduction = "pca", dims = 1:15)
sub <- FindClusters(sub, resolution = 0.2, algorithm = 2)
sub <- RunUMAP(sub, reduction = "pca", dims = 1:15, reduction.name = "umap")
DimPlot(sub,reduction = "umap",raster=FALSE)
sub$ALRA_clusters <- Idents(sub)
p <- DimPlot(sub, reduction = "umap", label = TRUE)
ggsave("/scratch/khandl/technical/figures/Dropouts/DimPlot_alra_clusters.svg", width = 10, height = 8, plot = p)

DimPlot(sub, group.by = "eos_subtypes")
markers <-  c("Mki67","Tuba1b","Epx","Prg2","Ear1","Ear2","Alox15","Aldh2","S100a9","S100a6","S100a10","Il5","Retnla","Ccl9","Il1rl1","Cd24a","Mmp9",
              "Icosl","Il4","Tgfb1","Pirb","Rara","Cd80","Cd274","Ptgs2","Il1rn","Tnf","Siglecf")
p <- DotPlot(sub, features = markers,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/Dropouts/RNA_alra_dotPlot_markers.svg", width = 10, height = 8, plot = p)

p <- DimPlot(sub, group.by = "eos_subtypes", reduction = "umap", label = TRUE, cols = c("#10A069","#E81818" ,"#878283","#E9EF11", "#26DFED","#E88A1A")) 
ggsave("/scratch/khandl/technical/figures/Dropouts/DimPlot_alra_eos_subtypes.svg", width = 10, height = 8, plot = p)

p <- FeaturePlot(sub, features = c("Epx","Ear1","Cd80","Cd274"), ncol = 2, cols = c('lightgrey','darkred'))
ggsave("/scratch/khandl/technical/figures/Dropouts/Eos_marker_alra.svg", width = 15, height = 20, plot = p)

### compare the clustering overlap with mclust
table(sub$RNA_clusters, sub$ALRA_clusters)

# Adjusted Rand Index (ARI) is a statistical measure that is used in cluster analysis to evaluate the similarity between two different clusterings 
ari <- adjustedRandIndex(sub$RNA_clusters, sub$ALRA_clusters)
ari

