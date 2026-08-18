########### This code analyses marker gene expression in eosinophils, neutrophils, mast cells and basophils ##########
### Datasets used: GSE256088, GSE175930, E-MTAB-14010, GSM7919060, GSE276583, GSE216189, GSE282765, GSE182001

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Mm
obj1 <- readRDS(file.path(seurat_objects_dir,"Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj2 <- readRDS(file.path(seurat_objects_dir,"Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj3 <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj4 <- readRDS(file.path(seurat_objects_dir,"Mm_bm_GSM7919060_anno.rds"))
obj5 <- readRDS(file.path(seurat_objects_dir,"Mm_liver_GSE216189.rds"))

obj <- merge(obj1, c(obj2,obj3,obj4,obj5))
obj <- JoinLayers(obj)

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B", "Basophils", "DCs", "EoP","Eosinophils","GMPs","HSCs", "Kupffer", "Macrophages","Mast",
                              "Monocytes","Neutrophils","PCs", "ProMono", "ProNeutro","T","TAMs"))

### Eos, mast, mono, baso, neutrophil markers 
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

Idents(obj) <- "annotation"
sub <- subset(obj, idents = c("Mast","Eosinophils","Neutrophils","Basophils"))

markers_oi <- c("Cpa3","Fcer1a","Ms4a2","Gata2", # Mast and Basophils
                "Mcpt8","Cd200r3","Prss34", # Basophils
                "Tpsab1", "Tpsb2","Kit", # Mast cells
                "Ccr3","F5","Syne1","Dach1","Ccl6","Adgre1","Arhgap18","Cyp4f18","Adam8", # Eosinophils
                "S100a8","S100a9","Csf3r", "Ngp","Ly6g" # Neutrophils
)
p <- DotPlot(sub, features = markers_oi,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_discrete(limits =c("Mast","Basophils","Eosinophils", "Neutrophils"))
ggsave(file.path(eos_marker_plots_dir, "Mm_markers1.svg"), width = 15, height = 6, plot = p)

markers_oi <- c("Ngp","Mmp9","Camp","Lcn2","S100a9","S100a8","Chil3","Lyz2"
)
sub$anno_cond <- paste0(sub$annotation,"_",sub$condition)
Idents(sub) <- "anno_cond"
p <- DotPlot(sub, features = markers_oi,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
  scale_y_discrete(limits =c("Mast","Basophils","Eosinophils", "Neutrophils"))
ggsave(file.path(eos_marker_plots_dir, "Mm_markers2.svg"), width = 15, height = 6, plot = p)

##### Hs
obj1 <- readRDS(file.path(seurat_objects_dir,"Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj2 <- readRDS(file.path(seurat_objects_dir,"Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds"))
obj3 <- readRDS(file.path(seurat_objects_dir,"Hs_blood_GSE256088_et_al_anno.rds"))
obj4 <- readRDS(file.path(seurat_objects_dir,"Hs_esophagus_duodenum_EoE_GSE175930_anno.rds"))
obj5 <- readRDS(file.path(seurat_objects_dir,"Hs_blood_E-MTAB-14010_anno.rds"))
obj6 <- readRDS(file.path(seurat_objects_dir,"Hs_blood_GSE276583_anno.rds"))

obj <- merge(obj1, c(obj2,obj3,obj4,obj5,obj6))
obj <- JoinLayers(obj)

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B","Basophils","DCs","Eosinophils","Macrophages","Mast","Monocytes","Neutrophils","PCs","ProNeutro", "T","TAMs"))

### Eos, mast, mono, baso, neutrophil markers 
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

sub2 <- subset(obj, idents = c("Mast","Eosinophils","Neutrophils","Basophils"))

markers_oi <- c("CPA3","MS4A2","FCER1A","HDC", # Basophils and Mast cells 
                "ENPP3","CD200R1","CLEC12A", # Basophils 
                "TPSAB1","KIT","TPSB2", #Mast cells 
                "CLC", "SMPD3","SYNE1","ALOX15","PRSS33", # Eosinophils 
                "CSF3R","FCGR3B","CXCR2", "MME","CXCL8" # Neutrophils
)
p <- DotPlot(sub2, features = markers_oi,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  scale_y_discrete(limits =c("Mast","Basophils","Eosinophils", "Neutrophils")) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave(file.path(eos_marker_plots_dir, "Hs_markers.svg"), width = 15, height = 6, plot = p)
