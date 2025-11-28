########### This code integrates eosinophils from WT and Il5-tg mice, colon and BM ##########
### Datasets used: GSE182001, GSE282765

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.7.Functions_cell_type_prop.R")

##### Load Seurat objects
il5tg <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")
wt_bm <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
wt_colon <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")

### Extract BM and colon 
Idents(wt_bm) <- "condition"
wt_bm <- subset(wt_bm, idents =c("bm_wt","blood_wt") )
Idents(wt_colon) <- "condition"
wt_colon <- subset(wt_colon, idents =c("adult_colon_wt") )
Idents(il5tg) <- "condition"
il5tg <- subset(il5tg, idents =c("bm","colon","blood","spleen") )

### Extract Eosinophils and EoP
Idents(wt_bm) <- "annotation"
wt_bm <- subset(wt_bm, idents = c("Eosinophils","EoP"))
Idents(wt_colon) <- "annotation"
wt_colon <- subset(wt_colon, idents = c("Eosinophils"))
Idents(il5tg) <- "annotation"
il5tg <- subset(il5tg, idents = c("Eosinophils","EoP"))

obj <- merge(wt_bm, y= c(wt_colon,il5tg))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

saveRDS(obj, "/scratch/khandl/technical/seurat_objects/Il5tg_wt_eos.rds" )

current.cluster.ids <- c("adult_colon_wt","bm","bm_wt","colon")
new.cluster.ids <- c("wt","il5tg","wt","il5tg")
obj$genotype <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

current.cluster.ids <- c("adult_colon_wt","bm","bm_wt","colon")
new.cluster.ids <- c("colon","bm","bm","colon")
obj$tissue <- plyr::mapvalues(x = obj$condition, from = current.cluster.ids, to = new.cluster.ids)

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

obj <- FindNeighbors(obj, reduction = "pca", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.8, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:15, reduction.name = "umap")
DimPlot(obj,reduction = "umap",raster=FALSE, label= TRUE, label.size = 8)
DimPlot(obj,reduction = "umap",raster=FALSE, label= TRUE, label.size = 8,group.by = "condition")
DimPlot(obj,reduction = "umap",raster=FALSE, label= TRUE, label.size = 8,group.by = "genotype")

obj <- JoinLayers(obj)

#### FastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$genotype)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:8)
obj <- FindClusters(obj, resolution = 1, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:8, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 10)
DimPlot(obj,reduction = "umap.mnn",group.by = "genotype", label = TRUE, label.size = 10)
DimPlot(obj,reduction = "umap.mnn",group.by = "condition", label = TRUE, label.size = 10)

DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label= TRUE, label.size = 8, split.by = "condition")

## check quality of clusters
VlnPlot(obj, features = "nFeature_RNA")
# remove cluster 8, has very low nFeaturs, even lower than the colon eosinophils 

Idents(obj) <- "mnn.clusters"
obj <- subset(obj, idents = c(0:7,9:11))

p <- DimPlot(obj, group.by = "tissue", reduction = "umap.mnn", label = TRUE, split.by = "genotype", cols = c("#26DFED","#E81818" )) 
ggsave("/scratch/khandl/technical/figures/wt_il5tg/umap_split_genotype_group_tissue.svg", width = 10, height = 8, plot = p)

### DotPlot of marker genes
goi <- c("Mki67","Tuba1b", "Epx","Prg3", #progenitor
         "Ear1","Ear2","Alox15","Cebpe", # immature 
         #"Aldh2","S100a6","Retnla","Ccl9","Il1rl1","Cd24a", # circulating 
         "Mmp9","Icosl","Il4","Tgfb1","Pirb","Rara","Krt80","Tmem71","Lcp2","Sell","Il17ra", # basal 
         "Cd80","Cd274","Ptgs2","Il1rn","Il1b","Vegfa" # active 
)
Idents(obj) <- "mnn.clusters"
DotPlot(obj, features = goi,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

# Rename
current.cluster.ids <- c(0:7,9:11)
new.cluster.ids <- c("active","active","active","active" ,"active","precursor","basal","basal","immature","immature","precursor")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn", split.by = "tissue")
obj_save <- obj

### DEGs per annotation cluster 
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

obj$anno_cond <- paste0(obj$annotation, "_",obj$genotype)

obj$anno_cond <- factor(obj$anno_cond, levels = c("precursor_il5tg","precursor_wt","immature_il5tg","immature_wt", "basal_il5tg","basal_wt","active_il5tg","active_wt"))
Idents(obj) <- "anno_cond"
p <- DotPlot(obj, features = goi,dot.scale = 10, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/technical/figures/wt_il5tg/Dotplot_annotation_marker_per_genotype.svg", width = 10, height = 8, plot = p)

p <- DimPlot(obj, group.by = "annotation", reduction = "umap.mnn", label = TRUE, split.by = "genotype", cols = c("#E81818" ,"#26DFED","#10A069","#E88A1A" )) 
ggsave("/scratch/khandl/technical/figures/wt_il5tg/umap_split_genotype_group_annotation.svg", width = 10, height = 8, plot = p)

##### Subtype proportions per genotype
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/technical/figures/wt_il5tg/","condition")
df <- read.csv("/scratch/khandl/technical/figures/wt_il5tg/condition_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2,3,4,6)) 

#define order 
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c( "precursor","immature","basal","active" )))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#26DFED","#E88A1A", "#10A069","#E81818" )) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/technical/figures/wt_il5tg/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### QC between il5-tg and wt 
obj$tissue_cond <- paste0(obj$tissue, "_",obj$genotype)
p <- VlnPlot(obj, features= "nFeature_RNA", group.by = "tissue_cond", pt.size = 0, cols = c("#26DFED","#6DA0D5" ,"#E81818", "#6D0A16" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  +  
  scale_x_discrete(limits =c("bm_il5tg", "bm_wt","colon_il5tg","colon_wt"))
print(p)
ggsave("/scratch/khandl/technical/figures/wt_il5tg/nFeature.svg", width = 8, height = 8, plot = p)

p <- VlnPlot(obj, features= "percent.mt", group.by = "tissue_cond", pt.size = 0, cols = c("#26DFED","#6DA0D5" ,"#E81818", "#6D0A16" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  +  
  scale_x_discrete(limits =c("bm_il5tg", "bm_wt" ,"colon_il5tg","colon_wt"))
print(p)
ggsave("/scratch/khandl/technical/figures/wt_il5tg/percent.mt.svg", width = 8, height = 8, plot = p)

##### Save object 
saveRDS(obj, "/scratch/khandl/technical/seurat_objects/wt_il5tg_bm_colon_annotated.rds")

