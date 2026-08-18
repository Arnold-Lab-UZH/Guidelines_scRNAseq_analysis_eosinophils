########### This code extracts eosinophil specific marker genes across datasets  ##########
### Datasets used: GSE256088, GSE175930, E-MTAB-14010, GSM7919060, GSE276583, GSE216189, GSE282765, GSE182001

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load data and calculate DEGs for eosinophils and EoP
### Hs
## NAT/tumor GSE282765 
obj <- readRDS(file.path(seurat_objects_dir,"Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# NAT 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Hs_NAT_GSE282765_MAST.csv"))

# tumor 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tumor")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Hs_tumor_GSE282765_MAST.csv"))

## PB GSE282765 
obj <- readRDS(file.path(seurat_objects_dir,"Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# PB   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE282765_MAST.csv"))

## PB  GSE256088 
obj <- readRDS(file.path(seurat_objects_dir,"Hs_blood_GSE256088_et_al_anno.rds"))
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE256088_MAST.csv"))

## Esophagus and duodenum GSE175930 
obj <- readRDS(file.path(seurat_objects_dir,"Hs_esophagus_duodenum_EoE_GSE175930_anno.rds"))

# Esophagus 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = c("esophagus"))
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Esophagus_GSE175930_MAST.csv"))

# Duodenum 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = c("duodenum"))
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_duodenum_GSE175930_MAST.csv"))

## PB 
obj <- readRDS(file.path(seurat_objects_dir,"Hs_blood_E-MTAB-14010_anno.rds"))
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_PB_E-MTAB-14010_MAST.csv"))

## PB GSE276583
obj <- readRDS(file.path(seurat_objects_dir,"Hs_blood_GSE276583_anno.rds"))
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE276583_MAST.csv"))

## Hs 
df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Hs_NAT_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Hs_tumor_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE256088_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Esophagus_GSE175930_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_duodenum_GSE175930_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:12, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_PB_E-MTAB-14010_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:12, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_PB_GSE276583_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

### Mm 
## tumor, colon GSE282765
obj <- readRDS(file.path(seurat_objects_dir,"Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# colon, NAT  
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_colon_GSE282765_MAST.csv"))

# tumor, disseminated  
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tumor")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars ="nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_tumor_GSE282765_MAST.csv"))

## BM and blood GSE282765 
obj <- readRDS(file.path(seurat_objects_dir,"Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# blood   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_blood_GSE282765_MAST.csv"))

# bm   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_bm_GSE282765_MAST.csv"))

## il5tg  GSE182001 
obj <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# bm   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_bm_GSE182001_MAST.csv"))

# blood   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_blood_GSE182001_MAST.csv"))

# spleen   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "spleen")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_spleen_GSE182001_MAST.csv"))

# stomach   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "stomach")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars ="nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_stomach_GSE182001_MAST.csv"))

# duodenum   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "small_int")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_duodenum_GSE182001_MAST.csv"))

# colon   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_colon_GSE182001_MAST.csv"))

## BM GSM7919060
obj <- readRDS(file.path(seurat_objects_dir,"Mm_bm_GSM7919060_anno.rds"))
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_bm_GSM7919060_MAST.csv"))

## liver GSE216189
obj <- readRDS(file.path(seurat_objects_dir,"Mm_liver_GSE216189.rds"))
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data",test.use = "MAST", latent.vars = "nFeature_RNA")
write.csv(markers, file.path(eos_marker_tables_dir, "Eos_marker_Mm_Liver_GSE216189_MAST.csv"))

## Mm
df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_colon_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:11, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_tumor_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_blood_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:11, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_bm_GSE282765_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "EoP",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_bm_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "EoP",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_blood_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_spleen_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_stomach_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:11, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_duodenum_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:11, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_il5tg_colon_GSE182001_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_bm_GSM7919060_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

df <- read.csv(file.path(eos_marker_tables_dir, "Eos_marker_Mm_Liver_GSE216189_MAST.csv"))
marker <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
marker[order(-marker$avg_log2FC), ][1:10, ]

