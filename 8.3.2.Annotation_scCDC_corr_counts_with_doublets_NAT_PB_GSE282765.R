########## This code uses scCDC for ambient RNA quantification  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.5.Functions_ambient_RNA_plotting.R"))

##### Load R object
obj <- readRDS(file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_annotated.rds"))

##### Run for each expriment/cartridge separate  
### Exp 1 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp1")

## Define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp1 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 2
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp2")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp2 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 3 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp3")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp3 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 4
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp4")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp4 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 5 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp5")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp5 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 7 = patient 6
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp7")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp6 = ContaminationCorrection(sub,rownames(GCGs))

### Exp 8 = patient 7
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp8")
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = file.path(ambient_rna_plots_dir, "GCGs/"))
rownames(GCGs)
exp7 = ContaminationCorrection(sub,rownames(GCGs))

##### Merge Seurat objects 
### overwrite the RNA assay with the corrected one 
exp1@assays$RNA <- exp1@assays$Corrected
exp1@assays$mnn.reconstructed <- NULL
exp2@assays$RNA <- exp2@assays$Corrected
exp2@assays$mnn.reconstructed <- NULL
exp3@assays$RNA <- exp3@assays$Corrected
exp3@assays$mnn.reconstructed <- NULL
exp4@assays$RNA <- exp4@assays$Corrected
exp4@assays$mnn.reconstructed <- NULL
exp5@assays$RNA <- exp5@assays$Corrected
exp5@assays$mnn.reconstructed <- NULL
exp6@assays$RNA <- exp6@assays$Corrected
exp6@assays$mnn.reconstructed <- NULL
exp7@assays$RNA <- exp7@assays$Corrected
exp7@assays$mnn.reconstructed <- NULL

merged <- merge(exp1, y = c(exp2, exp3,exp4,exp5,exp6,exp7))
merged <- JoinLayers(merged)

##### Clustering 
### Pre-processing 
DefaultAssay(merged) <- "RNA"
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$experiment)
merged <- NormalizeData(merged,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
merged <- RunPCA(merged, features = VariableFeatures(object =merged), npcs = 20, verbose = FALSE)

### FastMNN integration 
merged <- IntegrateLayers(object = merged, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                          verbose = FALSE)
ElbowPlot(merged)         
merged <- FindNeighbors(merged, reduction = "integrated.mnn", dims = 1:15)
merged <- FindClusters(merged, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
merged <- RunUMAP(merged, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(merged,reduction = "umap.mnn",group.by = "mnn.clusters", label = TRUE, label.size = 8)
merged <- JoinLayers(merged)

##### Cluster annotation 
### DEGs per cluster 
Idents(merged) <- "mnn.clusters"
markers <- FindAllMarkers(object = merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### nFeature and percent.mito per cluster to exclude low quality clusters 
VlnPlot(merged, features = "nFeature_RNA", pt.size = 0)
VlnPlot(merged, features = "percent.mt", pt.size = 0)

### Marker gene expression 
DotPlot(merged, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" # endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:25)
new.cluster.ids <- c("Neutrophils","T","PCs","Eosinophils", "Mono_DCs_Mac", "lowQ",
                     "Mono_DCs_Mac","lowQ","PCs", "Epithelial","Neutrophils",
                     "Endothelial","Mast","Fibroblasts","PCs" ,"B",
                     "PCs","T","Mono_DCs_Mac","Undefined","Mixed","Fibroblasts",
                     "Mast","Mixed","Mixed","Mixed")
merged$annotation <- plyr::mapvalues(x = merged$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(merged, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

Idents(merged) <- "annotation"
DotPlot(merged, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
                                 "KLK10","CD200R1","CLEC12A", # Basophils 
                                 "CPA3","MS4A2","FCER1A","ENPP3","HDC", # Basophils and Mast cells 
                                 "CLEC9A","FLT3","XCR1" ,# DCs
                                 "CD300E","EREG","VCAN", # Monocytes
                                 "F13A1", "FOLR2","C1QC" ,# Macrophages 
                                 "SPP1","MMP12","FN1", # TAMs
                                 "CLC","ADGRE1","SYNE1", # Eosinophils 
                                 "ICOS","CD8A", "CD3E","TRAC", # T cells 
                                 "CD19", "CR2","MS4A1", # B cells 
                                 "IGHG2","IGKC", # PCs
                                 "FCGR3B","S100A8", # Neutrophils 
                                 "EPCAM","TFF3","MUC2", # Epithelial 
                                 "COL1A2", # Fibroblasts/
                                 "PECAM1" # endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

### Subcluster Mono_DCs_Mac
Idents(merged) <- "annotation"
subCl <- FindSubCluster(merged,cluster = "Mono_DCs_Mac",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mono_DCs_Mac_0","Mono_DCs_Mac_1","Mono_DCs_Mac_2","Mono_DCs_Mac_3","Mono_DCs_Mac_4"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
FeaturePlot(sub_celltype, features = "CD1C")
FeaturePlot(sub_celltype, features = "VCAN")
# 3 and 4 = mixed, 2 = TAMs, 0 = Macs; 1 = mono_DC

# Rename
current.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Mast","Mixed",
                         "Mono_DCs_Mac_0","Mono_DCs_Mac_1","Mono_DCs_Mac_2","Mono_DCs_Mac_3","Mono_DCs_Mac_4",
                         "Neutrophils","PCs","T","Undefined")
new.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Mast","Mixed",
                     "Macrophages","Mono_DCs","TAMs","Mixed","Mixed",
                     "Neutrophils","PCs","T","Undefined")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster Mono_DCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mono_DCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mono_DCs_0","Mono_DCs_1","Mono_DCs_2"))
DimPlot(sub_celltype, reduction = "umap.mnn")

Idents(sub_celltype) <- "sub.cluster"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))
VlnPlot(sub_celltype, features = "nFeature_RNA", pt.size = 0)
FeaturePlot(sub_celltype, features = "CD1C")
FeaturePlot(sub_celltype, features = "VCAN")
# 0: DCs; 1, 2 = Monocytes

# Rename
current.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macrophages", "Mast","Mixed",
                         "Mono_DCs_0","Mono_DCs_1","Mono_DCs_2",
                         "Neutrophils","PCs","T","TAMs", "Undefined")
new.cluster.ids <- c("B","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages", "Mast","Mixed",
                     "DCs","Monocytes","Monocytes",
                     "Neutrophils","PCs","T","TAMs", "Undefined")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

##### Check annotation
Idents(subCl) <- "annotation"

DotPlot(subCl, features = unique(c("TPSAB1","KIT","TPSB2", # Mast cells 
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
                                   "FCGR3B","S100A8", # Neutrophils 
                                   "EPCAM","TFF3","MUC2", # Epithelial 
                                   "COL1A2", # Fibroblasts
                                   "PECAM1" # Endothelial 
)) )+ theme(axis.text.x = element_text(angle = 90)) 

##### Save object 
saveRDS(subCl, file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated.rds"))

##### control ambient RNA correction 
obj_scCDC <- readRDS(file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated.rds"))
obj_RNA <- readRDS(file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_annotated.rds"))

Idents(obj_scCDC) <- "annotation"
obj_scCDC <- subset(obj_scCDC, idents = c("PCs","Eosinophils"))
Idents(obj_RNA) <- "annotation"
obj_RNA <- subset(obj_RNA, idents = c("PCs","Eosinophils"))

### RNA assay 
DefaultAssay(obj_RNA) <- "RNA"
sce <- as.SingleCellExperiment(obj_RNA)

markers <- list(PC_Markers = c("IGKC","JCHAIN","IGLC2"),
                Eos_Marker = c("CLC","ALOX15","SYNE1"))

colData(sce)$cluster <- sce$annotation
df_RNA <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")

### scdc
DefaultAssay(obj_scCDC) <- "RNA"
obj_scCDC <- JoinLayers(obj_scCDC, assay = "RNA")
sce <- as.SingleCellExperiment(obj_scCDC, assay = "RNA")

markers <- list(PC_Markers = c("IGKC","JCHAIN","IGLC2"),
                Eos_Marker = c("CLC","ALOX15","SYNE1"))

colData(sce)$cluster <- sce$annotation
df_scdc <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")
df_scdc$assay <- "scdc"

df <- rbind(df_RNA, df_scdc)

## Plot marker labels of PC markers 
df1 <- df[df$markerLabels %in% "PC_Markers",]

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","scdc"="#13A2E2"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(doublet_plots_dir, "PCmarkers_after_sc.svg"), width = 10, height = 8, plot = p)

## Plot marker labels of Eosinophil markers 
df1 <- df[df$markerLabels %in% "Eos_Marker",]

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","scdc"="#13A2E2"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +  ylim(0,100)
ggsave(file.path(doublet_plots_dir, "EosMarker_after_sc.svg"), width = 10, height = 8, plot = p)


