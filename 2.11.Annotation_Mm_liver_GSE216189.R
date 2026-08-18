########## This code does sample integration, pre-processing, clustering and annotation of Mm liver from healthy mice and liver metastasis (AKPS) from GSE216189  ##########

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir,"1.2.Functions_Seurat_integration.R"))

##### Loading of Mouse gene/ensemble IDs to convert ensemble IDs to gene IDs 
ensembl<-useEnsembl(biomart="ensembl")
list<-listDatasets(ensembl)
mart <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",version=100)
attributes<-listAttributes(mart)
gene_ids_mouse <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"), mart = mart)
write.csv(gene_ids_mouse, file.path(annotation_tables_dir, "gene_ids_ensemble_ids_mouse.csv"))

##### Seurat object generation 
liver_mets1 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM7026580_Mets_sc1.dgecounts.rds"), 
                                          project = "liverMets",3,200, condition = "Liver_mets1",batch = "batch1",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets2 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM7026581_Mets_sc2.dgecounts.rds"), 
                                          project = "liverMets", 3,200,  condition = "Liver_mets2",batch = "batch1",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets3 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661394_WT.dgecounts.rds"), 
                                          project = "liverMets",3,200,   condition = "Liver_healthy1",batch = "batch2",tissue_oi = "liver",phenotype_oi = "healthy")

liver_mets4 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661376_MetsM1.dgecounts.rds"), 
                                          project = "liverMets",3,200,    condition = "Liver_mets3",batch = "batch3",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets5 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661378_MetsM3.dgecounts.rds"), 
                                          project = "liverMets",3,200,  condition = "Liver_mets4",batch = "batch3",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets6 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661380_MetsM4.dgecounts.rds"), 
                                          project = "liverMets",3,200,  condition = "Liver_mets5",batch = "batch3",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets7 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661382_Mets4M1.dgecounts.rds"), 
                                          project = "liverMets",3,200, condition = "Liver_mets6",batch = "batch4",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets8 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661384_Mets4M2.dgecounts.rds"), 
                                          project = "liverMets", 3,200, condition = "Liver_mets7",batch = "batch4",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets9 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661388_Mets6M1.dgecounts.rds"), 
                                          project = "liverMets",3,200,  condition = "Liver_mets8",batch = "batch5",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets10 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661390_Mets6M2.dgecounts.rds"), 
                                           project = "liverMets",3,200,  condition = "Liver_mets9",batch = "batch5",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets11 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661392_Mets6M3.dgecounts.rds"), 
                                           project = "liverMets",3,200,  condition = "Liver_mets10",batch = "batch5",tissue_oi = "liver",phenotype_oi = "metastasis")

liver_mets12 <- create_seurat_zUMI_outputs(file.path(raw_data_GSE216189_dir,"GSM6661386_MetsM5.dgecounts.rds"), 
                                           project = "liverMets",3,200,  condition = "Liver_mets11",batch = "batch6",tissue_oi = "liver",phenotype_oi = "metastasis")

### Merge samples
obj1 <- merge(liver_mets1, y = c(liver_mets2,liver_mets3,liver_mets4,liver_mets5,liver_mets6,liver_mets7,liver_mets8,liver_mets9,
                                 liver_mets10,liver_mets11,liver_mets12),
              add.cell.ids = c("1","2","3","4","5","6","7","8","9","10","11","12"))
obj1 <- JoinLayers(obj1)

### Add mitochondrial percentage per cell 
obj1$percent.mt <- PercentageFeatureSet(obj1, pattern = "^mt-")

### Apply mitochondrial cutoff 
obj1 <- subset(obj1, subset = nFeature_RNA < 8000 & percent.mt < 25)

##### Clustering 
### Pre-processing 
obj <- obj1
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

### FastMNN integration 
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
ElbowPlot(obj)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn",return.model = TRUE)
DimPlot(obj,reduction = "umap.mnn", label = TRUE, group.by = "mnn.clusters", label.size = 5)
DimPlot(obj,reduction = "umap.mnn", label = TRUE, group.by = "batch")
obj <- JoinLayers(obj)

##### Cluster annotation 
### DEGs per cluster 
DimPlot(obj,reduction = "umap.mnn", group.by = "mnn.clusters",label = TRUE)
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

### nFeature and percent.mito per cluster to exclude low quality clusters 
Idents(obj) <- "mnn.clusters"
VlnPlot(obj, features = "nFeature_RNA",pt.size = 0)
VlnPlot(obj, features = "percent.mt",pt.size = 0)

### Marker gene expression 
markers  <- c("Ighm", "Scd1","Cd19", # B cells
              "Jchain","Igha","Igkc", # PCs
              
              "Trac","Cd3d" ,"Icos","Cd8a","Klra8", # T
              
              "Flt3","Xcer1", "Clec9a", # cDCs
              "S100a4", "Itgax","Ly6c2","Vcan", # Monocytes
              "C1qc","Dnase1l3","Folr2","Lyve1", # Macrophages 
              "Spp1","Fn1","Mmp12","Arg1", # TAMs
              "Clec4f","Vsig4", # Kupffer cells
              
              "Hba-a2", # Red blood cells
              
              "Csf3r", "S100a8","S100a9", # eutrophils
              "Tpsab1", "Tpsb2","Kit", # Mast cells
              "Mcpt8","Cd200r3","Clec12a", # Basophils
              "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
              "Siglecf","Ccr3","Alox15","F5","Syne1", # Eosinophils 
              
              "Glul","Cyp2e1", "Ass1","Alb", # Hepatocytes

               "Pecam1","Dll4","Galnt15",  "Plpp1", # (LECs) Endothelial 
              "Lrat","Reln", # Stellate cells
              "Carmn","Nr1h5",  # Stromal cells liver
              "Svep1","Ncam1", # Fibroblasts liver

              "Sprr2a2","Pglyrp1","Gpx2","Epcam" # Epithelial 
)              

DotPlot(obj, features = unique(markers),dot.scale = 6, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

### Rename clusters
current.cluster.ids <- c(0:21)
new.cluster.ids <- c("Endothelial","Kupffer","B","T","Kupffer","Monocytes","lowQ","Fibroblasts", "DCs","Endothelial",
                     "Hepatocytes","Neutrophils","Mast_Eos_Baso","Epithelial", "DCs","?","?","PCs","?","?","?","?")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj,reduction = "umap.mnn", label = TRUE, group.by = "annotation")

### Subcluster Endothleial 
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Endothelial",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Endothelial_0","Endothelial_1","Endothelial_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c(
  "Pecam1","Dll4", # LSECs 
  "Plpp1", # Pericentral LSECs
  "Galnt15" # Periportal LSECs
)              

DotPlot(sub_celltype, features = unique(markers),dot.scale = 6, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial_0","Endothelial_1","Endothelial_2","Epithelial", "Hepatocytes","Kupffer","lowQ","Mast_Eos_Baso",
                         "Monocytes","Neutrophils","PCs","Fibroblasts","T")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Endothelial","lowQ","Epithelial", "Hepatocytes","Kupffer","lowQ","Mast_Eos_Baso",
                     "Monocytes","Neutrophils","PCs","Fibroblasts","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster Kupffer 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Kupffer",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Kupffer_0","Kupffer_1","Kupffer_2","Kupffer_3","Kupffer_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c( 
  "Clec9a","Xcer1", "Tlr3", # cDC1
  "Cd209a","Flt3","Ccl22", # cDC2
  "Siglech","Cox6a2", #pDCs
  "S100a4", "Itgax","Ly6c2","Vcan", # Monocytes
  "C1qc","Dnase1l3","Folr2","Lyve1", # Macrophages 
  "Spp1","Fn1","Mmp12","Arg1", # TAMs
  "Clec4f","Vsig4","Vcam1", # Kupffer cells, Vcam1 high = periportal 
  
  "Pecam1","Dll4", # LSECs 
  "Plpp1", # Pericentral LSECs
  "Galnt15", # Periportal LSECs
  
  "Krt19", "Epcam","Sox9", # Cholangiocytes
  "Sprr2a2","Pglyrp1","Gpx2", # Metastatic cells
  
  "Lrat","Reln", # Stellate cells
  "Carmn","Nr1h5",  # Stromal cells liver
  "Svep1","Ncam1", # Fibroblasts liver
  
  "Epcam","Spp1", "Mmp7","Thbs1" # AKPS cancer cells tumor 
)              

DotPlot(sub_celltype, features = unique(markers),dot.scale = 6, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Epithelial", "Hepatocytes","Kupffer_0","Kupffer_1","Kupffer_2","Kupffer_3","Kupffer_4",
                         "lowQ","Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","Kupffer","Kupffer","Kupffer","lowQ",
                     "lowQ","Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster T cells
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_0","T_1","T_2","T_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c(
  "Trac","Cd3d" ,"Icos","Ctla4", # CD4 T cells
  "Ccl5","Cd8a","Ms4a4b", # CD8 T cells 
  "Klra8","Cma1" ,#N K cells
  "Il17rb","Hlf" # ILCs
  
)              

DotPlot(sub_celltype, features = unique(markers),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

FeaturePlot(sub_celltype, features = c("Cd3g","Cd8a","Icos","Ctla4"), reduction = "umap.mnn")

# Rename
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                         "Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T_0","T_1","T_2","T_3")
new.cluster.ids <- c("?", "B","DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                     "Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T","T","T","lowQ")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster DCs 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "DCs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "DCs_0","DCs_1","DCs_2","DCs_3","DCs_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c(  
  "Clec9a","Xcer1", "Tlr3", # cDC1
  "Cd209a","Flt3","Ccl22", # cDC2
  "Siglech","Cox6a2", # pDCs
  "S100a4", "Itgax","Ly6c2","Vcan", # Monocytes
  "C1qc","Dnase1l3","Folr2","Lyve1", # Macrophages 
  "Spp1","Fn1","Mmp12","Arg1", # TAMs
  "Clec4f","Vsig4","Vcam1" # Kupffer cells, Vcam1 high = periportal 
  
)              

DotPlot(sub_celltype, features = unique(markers),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "B","DCs_0","DCs_1","DCs_2","DCs_3","DCs_4", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                         "Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T")
new.cluster.ids <- c("?", "B", "DCs","Monocytes","DCs","DCs","DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                      "Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster Fibroblasts 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Fibroblasts",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Fibroblasts_0","Fibroblasts_1","Fibroblasts_2","Fibroblasts_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c( "Csf3r", "S100a8","S100a9", # Neutrophils:
               "Fcer1a","Cpa3","Mcpt4","Tpsab1", # Mast 
               "Bpi","Ceacam8","Ltf","Camp","Ms4a3", # Basophils 
               "Siglecf","Ccr3","Alox15","F5","Syne1", #Eosinophils 
               
               "Pecam1","Dll4", #LSECs 
               "Plpp1", # Pericentral LSECs
               "Galnt15", # Periportal LSECs
               
               "Krt19", "Epcam","Sox9", #C holangiocytes
               "Sprr2a2","Pglyrp1","Gpx2", #Metastatic cells
               
               "Lrat","Reln", # Stellate cells
               "Carmn","Nr1h5",  # Stromal cells liver
               "Svep1","Ncam1", # Fibroblasts liver
               
               "Epcam","Spp1", "Mmp7","Thbs1" # AKPS cancer cells tumor 
)              

DotPlot(sub_celltype, features = unique(markers),dot.scale = 6, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "B", "DCs",  "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                         "Mast_Eos_Baso","Monocytes","Neutrophils","PCs","Fibroblasts_0","Fibroblasts_1","Fibroblasts_2","Fibroblasts_3","T")
new.cluster.ids <- c("?", "B", "DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ_cytoplasma_lost",
                     "Mast_Eos_Baso","Monocytes","Neutrophils","PCs",
                     "lowQ","lowQ","lowQ","Stellate", "T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster monocytes 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Monocytes",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Monocytes_0","Monocytes_1","Monocytes_2","Monocytes_3","Monocytes_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")
FeaturePlot(sub_celltype, reduction = "umap.mnn",features = c("C1qc","Lyz2","Spp1","Fn1","Ly6c1"))

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

# Rename
current.cluster.ids <- c("?", "B", "DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                         "Mast_Eos_Baso",
                         "Monocytes_0","Monocytes_1","Monocytes_2","Monocytes_3","Monocytes_4"
                         ,"Neutrophils","PCs","Stellate", "T")
new.cluster.ids <- c("?", "B", "DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ",
                      "Mast_Eos_Baso", "Macrophages","Monocytes","Monocytes","lowQ","lowQ","Neutrophils","PCs","Stellate", "T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Subcluster Mast_Eos_Baso
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mast_Eos_Baso",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Mast_Eos_Baso_0","Mast_Eos_Baso_1","Mast_Eos_Baso_2","Mast_Eos_Baso_3","Mast_Eos_Baso_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

markers  <- c( "Csf3r", "S100a8","S100a9", # Neutrophils
              "Tpsab1", "Tpsb2","Kit", # Mast cells
              "Mcpt8","Cd200r3","Clec12a", # Basophils
              "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
              "Siglecf","Ccr3","Alox15","F5","Syne1"# Eosinophils 
)              

DotPlot(sub_celltype, features = markers, scale = FALSE) + theme(axis.text.x = element_text(angle = 90)) 


FeaturePlot(sub_celltype, reduction = "umap.mnn",features = c("Siglecf","Cpa3","Alox15","F5"))

markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
# 1 = Eos, 3 = Baso, rest is ? 

# Rename
current.cluster.ids <- c("?", "B", "DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ","lowQ_cytoplasma_lost","Macrophages",
                         "Mast_Eos_Baso_0","Mast_Eos_Baso_1","Mast_Eos_Baso_2","Mast_Eos_Baso_3","Mast_Eos_Baso_4","Monocytes",
                         "Neutrophils","PCs","Stellate", "T")
new.cluster.ids <- c("?", "B", "DCs", "Endothelial","Epithelial", "Hepatocytes", "Kupffer","lowQ","lowQ","Macrophages",
                     "?","Eosinophils","?","Basophils","?","Monocytes",
                     "Neutrophils","PCs","Stellate", "T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,reduction = "umap.mnn")

### Check annotation 
markers  <- c("Ighm", "Scd1","Cd19", # B cells
              "Jchain","Igha","Igkc", # PCs
              
              "Trac","Cd3d" ,"Icos","Cd8a","Klra8", # T
              
              "Flt3","Xcer1", "Clec9a", # cDCs
              "S100a4", "Itgax","Ly6c2","Vcan", # Monocytes
              "C1qc","Dnase1l3","Folr2","Lyve1", # Macrophages 
              "Spp1","Fn1","Mmp12","Arg1", # TAMs
              "Clec4f","Vsig4", # Kupffer cells
              
              "Hba-a2", # Red blood cells
              
              "Csf3r", "S100a8","S100a9", # eutrophils
              "Tpsab1", "Tpsb2","Kit", # Mast cells
              "Mcpt8","Cd200r3","Clec12a", # Basophils
              "Cpa3","Fcer1a","Ms4a2","Hdc", # Mast and Basophils
              "Siglecf","Ccr3","Alox15","F5","Syne1", # Eosinophils 
              
              "Glul","Cyp2e1", "Ass1","Alb", # Hepatocytes
              
              "Pecam1","Dll4","Galnt15",  "Plpp1", # (LECs) Endothelial 
              "Lrat","Reln", # Stellate cells
              "Carmn","Nr1h5",  # Stromal cells liver
              "Svep1","Ncam1", # Fibroblasts liver
              
              "Sprr2a2","Pglyrp1","Gpx2","Epcam" # Epithelial 
)              
DotPlot(subCl, features = unique(markers),dot.scale = 6, scale = TRUE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

##### Add conditions to metadata 
subCl$technology <- "BD"
subCl$species <- "Mm"
subCl$cell_enrichment <- "CD45pos"
subCl$tissue <- "liver"

##### Save object 
saveRDS(subCl,file=file.path(seurat_objects_dir,"Mm_liver_GSE216189.rds"))

