############ This code compares outputs of empty droplet characterisation from forced cell determination, EmptyDrops and MALAT1 ##########

# Hs tumor and NAT GSE282765
# take forced cell determination and intronic + exonic reads

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load objects 
obj_forced_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
obj_automatic_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_automatic_cell_determination_with_intronic_reads_annotated.rds")
obj_emptyDrops_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_emptyDrops_determination_with_intronic_reads_annotated.rds")
obj_Malat1_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_MALAT1_determination_with_intronic_reads_annotated.rds")
obj_with_low200_high25 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_lower200features_larger25mito.rds")

## extract eosinophils
Idents(obj_forced_all) <- "annotation"
obj_forced <- subset(obj_forced_all, idents = "Eosinophils")
Idents(obj_automatic_all) <- "annotation"
obj_automatic <- subset(obj_automatic_all, idents = "Eosinophils")
Idents(obj_emptyDrops_all) <- "annotation"
obj_emptyDrops <- subset(obj_emptyDrops_all, idents = "Eosinophils")
Idents(obj_Malat1_all) <- "annotation"
obj_Malat1 <- subset(obj_Malat1_all, idents = "Eosinophils")
Tcells <- subset(obj_forced_all, idents = "T")

obj_with_low200_high25$cell_determination <- obj_with_low200_high25$annotation
Tcells$cell_determination <- "T"
##### for forced, Malat1 and emptyDrops subset cell IDs that are not present in automatic 
automatic_cell_ids <- rownames(obj_automatic@meta.data)
forced_cell_ids <- rownames(obj_forced@meta.data)
emptyDrops_cell_ids <- rownames(obj_emptyDrops@meta.data)
Malat1_cell_ids <- rownames(obj_Malat1@meta.data)

forced_additional <- forced_cell_ids[!forced_cell_ids %in% intersect(automatic_cell_ids, forced_cell_ids)]
emptyDrops_additional <- emptyDrops_cell_ids[!emptyDrops_cell_ids %in% intersect(automatic_cell_ids, emptyDrops_cell_ids)]
Malat1_additional <- Malat1_cell_ids[!Malat1_cell_ids %in% intersect(automatic_cell_ids, Malat1_cell_ids)]

obj_forced <- subset(obj_forced, cells = forced_additional)
obj_emptyDrops <- subset(obj_emptyDrops, cells = emptyDrops_additional)
obj_Malat1 <- subset(obj_Malat1, cells = Malat1_additional)

### merge all together 
obj <- merge(obj_with_low200_high25, y = c(obj_automatic,obj_forced,obj_emptyDrops,obj_Malat1,Tcells), add.cell.ids = c("a","b","c","d","e","f"))
obj <- JoinLayers(obj)

##### marker gene expression 
markers <- c("CLC","CCR3","ALOX15","ADGRE1","DACH1", "SYNE1","SYNE2", "ADAM8", "MARCHF3","FFAR2","SSH2","PLAUR","CCL4L2","AREG","NFKB1","RELB","NCF2","CAMK1D", "OLR1","CD3E","ICOS","CD8A","CD4")

Idents(obj) <- "cell_determination"
p <- DotPlot(obj, features = markers,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits =c("features_lower200","mito_larger25","T","MALAT1" ,"EmptyDrops", "forced","automatic"))
ggsave("/scratch/khandl/technical/figures/summary_empty/DotPlot_markers_real_different_approaches.svg", width = 15, height = 6, plot = p)

p <- VlnPlot(obj, features= "nFeature_RNA", group.by = "cell_determination", pt.size = 0, cols = c("#7251A1","#6DA0D5" ,"#6D7272","#F4C914" ,"#23803F", "#020202","#5BC8D9" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  +  
  scale_x_discrete(limits =c("automatic", "forced","EmptyDrops","MALAT1" ,"T","mito_larger25","features_lower200"))
print(p)
ggsave("/scratch/khandl/technical/figures/summary_empty/nFeature.svg", width = 8, height = 8, plot = p)

p <- VlnPlot(obj, features= "percent.mt", group.by = "cell_determination", pt.size = 0, cols = c("#7251A1","#6DA0D5" ,"#6D7272","#F4C914" ,"#23803F", "#020202","#5BC8D9" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  + 
  scale_x_discrete(limits =c("automatic", "forced","EmptyDrops","MALAT1" ,"T","mito_larger25","features_lower200"))
print(p)
ggsave("/scratch/khandl/technical/figures/summary_empty/percent_mt.svg", width = 8, height = 8, plot = p)

