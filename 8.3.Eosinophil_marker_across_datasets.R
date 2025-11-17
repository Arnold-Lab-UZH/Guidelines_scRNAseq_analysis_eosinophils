########## This code identifies eosinophil specific markers ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load data and calculate DEGs for eosinophils and EoP
### Hs
## NAT/tumor GSE282765 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# NAT 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_GSE282765.csv")

# tumor 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tumor")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_tumor_GSE282765.csv")

## PB GSE282765 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# PB   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE282765.csv")

## PB  GSE256088 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_GSE256088_et_al_anno.rds")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE256088.csv")

## Esophagus and duodenum GSE175930 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_esophagus_duodenum_EoE_GSE175930_anno.rds")

# Esophagus 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = c("esophagus"))
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Esophagus_GSE175930.csv")

# Duodenum 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = c("duodenum"))
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_duodenum_GSE175930.csv")

## PB 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_E-MTAB-14010_anno.rds")

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_E-MTAB-14010.csv")

## PB GSE276583
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_GSE276583_anno.rds")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE276583.csv")

## Hs 
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_GSE282765.csv")
eos_marker_Hs_NAT_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_NAT_GSE282765 <- eos_marker_Hs_NAT_GSE282765[order(-eos_marker_Hs_NAT_GSE282765$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_tumor_GSE282765.csv")
eos_marker_Hs_CRC_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_CRC_GSE282765 <- eos_marker_Hs_CRC_GSE282765[order(-eos_marker_Hs_CRC_GSE282765$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE282765.csv")
eos_marker_Hs_PB_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_PB_GSE282765 <- eos_marker_Hs_PB_GSE282765[order(-eos_marker_Hs_PB_GSE282765$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE256088.csv")
eos_marker_Hs_PB_GSE256088 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_PB_GSE256088 <- eos_marker_Hs_PB_GSE256088[order(-eos_marker_Hs_PB_GSE256088$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Esophagus_GSE175930.csv")
eos_marker_Hs_Esophagus_GSE175930 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_Esophagus_GSE175930 <- eos_marker_Hs_Esophagus_GSE175930[order(-eos_marker_Hs_Esophagus_GSE175930$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_duodenum_GSE175930.csv")
eos_marker_Hs_duodenum_GSE175930 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_duodenum_GSE175930 <- eos_marker_Hs_duodenum_GSE175930[order(-eos_marker_Hs_duodenum_GSE175930$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_E-MTAB-14010.csv")
eos_marker_Hs_PB_E_MTAB_14010 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_PB_E_MTAB_14010 <- eos_marker_Hs_PB_E_MTAB_14010[order(-eos_marker_Hs_PB_E_MTAB_14010$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_PB_GSE276583.csv")
eos_marker_Hs_PB_GSE276583 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Hs_PB_GSE276583 <- eos_marker_Hs_PB_GSE276583[order(-eos_marker_Hs_PB_GSE276583$avg_log2FC), ][1:200, ]

## Pb and tissue 
intersect_PB_and_tissue_eos <- Reduce(intersect, list(eos_marker_Hs_NAT_GSE282765$gene, eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_PB_GSE282765$gene,
                                                      eos_marker_Hs_PB_GSE256088$gene,
                                                      eos_marker_Hs_Esophagus_GSE175930$gene, 
                                                      eos_marker_Hs_duodenum_GSE175930$gene, 
                                                      eos_marker_Hs_PB_E_MTAB_14010$gene,
                                                      eos_marker_Hs_PB_GSE276583$gene))

### Tissue Eos 
intersect_tissue_eos <- Reduce(intersect, list(eos_marker_Hs_NAT_GSE282765$gene, eos_marker_Hs_CRC_GSE282765$gene,
                                               eos_marker_Hs_Esophagus_GSE175930$gene, 
                                               eos_marker_Hs_duodenum_GSE175930$gene))

x <- list("NAT_GSE282765" = eos_marker_Hs_NAT_GSE282765$gene, "CRC_GSE282765" =eos_marker_Hs_CRC_GSE282765$gene,
          "Esophagus_GSE175930" = eos_marker_Hs_Esophagus_GSE175930$gene, "duodenum_GSE175930" = eos_marker_Hs_duodenum_GSE175930$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#EABBC4")

# unique per condition
setdiff(eos_marker_Hs_NAT_GSE282765$gene,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_duodenum_GSE175930$gene)))
setdiff(eos_marker_Hs_CRC_GSE282765$gene,Reduce(union, list(eos_marker_Hs_NAT_GSE282765$gene,eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_duodenum_GSE175930$gene)))
setdiff(eos_marker_Hs_Esophagus_GSE175930$gene,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_NAT_GSE282765$gene, eos_marker_Hs_duodenum_GSE175930$gene)))
setdiff(eos_marker_Hs_duodenum_GSE175930$gene,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_NAT_GSE282765$gene, eos_marker_Hs_Esophagus_GSE175930$gene)))

# shared NAT tumor 
shared_NAT_tumor <- intersect(eos_marker_Hs_NAT_GSE282765$gene,eos_marker_Hs_CRC_GSE282765$gene)
setdiff(shared_NAT_tumor,Reduce(union, list(eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_duodenum_GSE175930$gene)))

# shared Esophagus and duodenum 
shared_eso_duo <- intersect(eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_duodenum_GSE175930$gene)
setdiff(shared_eso_duo,Reduce(union, list(eos_marker_Hs_NAT_GSE282765$gene,eos_marker_Hs_CRC_GSE282765$gene)))

# shared Esophagus and tumor 
shared_eso_tumor <- intersect(eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_CRC_GSE282765$gene)
setdiff(shared_eso_tumor,Reduce(union, list(eos_marker_Hs_NAT_GSE282765$gene,eos_marker_Hs_duodenum_GSE175930$gene)))

# shared esophagus NAT 
shared_eso_NAT <- intersect(eos_marker_Hs_Esophagus_GSE175930$gene, eos_marker_Hs_NAT_GSE282765$gene)
setdiff(shared_eso_NAT,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_duodenum_GSE175930$gene)))

# shared duodenum NAT 
shared_duo_NAT <- intersect(eos_marker_Hs_duodenum_GSE175930$gene, eos_marker_Hs_NAT_GSE282765$gene)
setdiff(shared_duo_NAT,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_Esophagus_GSE175930$gene)))

# shared duodenum tumor
intersect_duo_tum <- intersect(eos_marker_Hs_duodenum_GSE175930$gene, eos_marker_Hs_CRC_GSE282765$gene)
setdiff(intersect_duo_tum,Reduce(union, list(eos_marker_Hs_NAT_GSE282765$gene,eos_marker_Hs_Esophagus_GSE175930$gene)))

# shared Esophagus, NAT, tumor 
shared_eos_NAT_tumor <- Reduce(intersect, list(eos_marker_Hs_Esophagus_GSE175930$gene,eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_NAT_GSE282765$gene))
setdiff(shared_eos_NAT_tumor,Reduce(union, list(eos_marker_Hs_duodenum_GSE175930$gene)))

# shared duodenum, esophagus, NAT 
shared_duo_eso_NAT <- Reduce(intersect, list(eos_marker_Hs_Esophagus_GSE175930$gene,eos_marker_Hs_duodenum_GSE175930$gene,eos_marker_Hs_NAT_GSE282765$gene))
setdiff(shared_duo_eso_NAT,Reduce(union, list(eos_marker_Hs_CRC_GSE282765$gene)))

# shared tumor, duodenum, NAT
shared_tum_duo_NAT <- Reduce(intersect, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_duodenum_GSE175930$gene,eos_marker_Hs_NAT_GSE282765$gene))
setdiff(shared_tum_duo_NAT,Reduce(union, list(eos_marker_Hs_Esophagus_GSE175930$gene)))

# shared tumor, esophagus, duodenum 
shared_tum_eos_duod <- Reduce(intersect, list(eos_marker_Hs_CRC_GSE282765$gene,eos_marker_Hs_duodenum_GSE175930$gene,eos_marker_Hs_Esophagus_GSE175930$gene))
setdiff(shared_tum_eos_duod,Reduce(union, list(eos_marker_Hs_NAT_GSE282765$gene)))

### Blood Eos 
intersect <- Reduce(intersect, list(eos_marker_Hs_PB_GSE282765$gene,
                                    eos_marker_Hs_PB_GSE256088$gene,
                                    eos_marker_Hs_PB_E_MTAB_14010$gene,
                                    eos_marker_Hs_PB_GSE276583$gene))

## venn diagram 
x <- list("PB_GSE282765" = eos_marker_Hs_PB_GSE282765$gene, "PB_GSE256088" =eos_marker_Hs_PB_GSE256088$gene,
          "PB_E_MTAB_14010" = eos_marker_Hs_PB_E_MTAB_14010$gene, "PB_GSE276583" = eos_marker_Hs_PB_GSE276583$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#B5DCB2")

# unique per condition
setdiff(eos_marker_Hs_PB_GSE282765$gene,Reduce(union, list(eos_marker_Hs_PB_GSE256088$gene,eos_marker_Hs_PB_E_MTAB_14010$gene, eos_marker_Hs_PB_GSE276583$gene)))
setdiff(eos_marker_Hs_PB_GSE256088$gene,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene,eos_marker_Hs_PB_E_MTAB_14010$gene, eos_marker_Hs_PB_GSE276583$gene)))
setdiff(eos_marker_Hs_PB_E_MTAB_14010$gene,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene,eos_marker_Hs_PB_GSE256088$gene, eos_marker_Hs_PB_GSE276583$gene)))
setdiff(eos_marker_Hs_PB_GSE276583$gene,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene,eos_marker_Hs_PB_GSE256088$gene, eos_marker_Hs_PB_E_MTAB_14010$gene)))

# shared 65 AND 88
intersect <- intersect(eos_marker_Hs_PB_GSE282765$gene,eos_marker_Hs_PB_GSE256088$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_E_MTAB_14010$gene, eos_marker_Hs_PB_GSE276583$gene)))

# shared 88 AND 10
intersect <- intersect(eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE256088$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene, eos_marker_Hs_PB_GSE276583$gene)))

# shared 10 AND 83
intersect <- intersect(eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE276583$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene, eos_marker_Hs_PB_GSE256088$gene)))

# shared 88 and 83 
intersect <- intersect(eos_marker_Hs_PB_GSE276583$gene,eos_marker_Hs_PB_GSE256088$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_E_MTAB_14010$gene, eos_marker_Hs_PB_GSE282765$gene)))

# shared 65 AND 83
intersect <- intersect(eos_marker_Hs_PB_GSE276583$gene,eos_marker_Hs_PB_GSE282765$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_E_MTAB_14010$gene, eos_marker_Hs_PB_GSE256088$gene)))

# shared 65 AND 10
intersect <- intersect(eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE282765$gene)
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE276583$gene, eos_marker_Hs_PB_GSE256088$gene)))

# shared 88, 10, 83
intersect <- Reduce(intersect, list(eos_marker_Hs_PB_GSE276583$gene,eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE256088$gene))
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE282765$gene)))

# shared 88, 10, 65
intersect <- Reduce(intersect, list(eos_marker_Hs_PB_GSE256088$gene,eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE282765$gene))
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE276583$gene)))

# shared 65,10,83
intersect <- Reduce(intersect, list(eos_marker_Hs_PB_GSE276583$gene,eos_marker_Hs_PB_E_MTAB_14010$gene,eos_marker_Hs_PB_GSE282765$gene))
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_GSE256088$gene)))

# shared 65,88,83
intersect <- Reduce(intersect, list(eos_marker_Hs_PB_GSE276583$gene,eos_marker_Hs_PB_GSE256088$gene,eos_marker_Hs_PB_GSE282765$gene))
setdiff(intersect,Reduce(union, list(eos_marker_Hs_PB_E_MTAB_14010$gene)))

### Mm 
## tumor, colon GSE282765
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# colon, NAT  
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_colon_GSE282765.csv")

# tumor, disseminated  
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tumor")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_tumor_GSE282765.csv")

## BM and blood GSE282765 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# blood   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_blood_GSE282765.csv")

# bm   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_bm_GSE282765.csv")

## il5tg  GSE182001 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

# bm   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_bm_GSE182001.csv")

# blood   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_blood_GSE182001.csv")

# spleen   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "spleen")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_spleen_GSE182001.csv")

# stomach   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "stomach")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_stomach_GSE182001.csv")

# duodenum   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "small_int")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_duodenum_GSE182001.csv")

# colon   
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")
sub <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_colon_GSE182001.csv")

## BM GSM7919060
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_bm_GSM7919060_anno.rds")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_bm_GSM7919060.csv")

## liver GSE216189
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Liver_GSE216189.rds")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_Liver_GSE216189.csv")

## Mm
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_colon_GSE282765.csv")
eos_marker_Mm_colon_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_colon_GSE282765 <- eos_marker_Mm_colon_GSE282765[order(-eos_marker_Mm_colon_GSE282765$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_tumor_GSE282765.csv")
eos_marker_Mm_tumor_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_tumor_GSE282765 <- eos_marker_Mm_tumor_GSE282765[order(-eos_marker_Mm_tumor_GSE282765$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_blood_GSE282765.csv")
eos_marker_Mm_blood_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_blood_GSE282765 <- eos_marker_Mm_blood_GSE282765[order(-eos_marker_Mm_blood_GSE282765$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_bm_GSE282765.csv")
eos_marker_Mm_bm_GSE282765 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "EoP",]
eos_marker_Mm_bm_GSE282765 <- eos_marker_Mm_bm_GSE282765[order(-eos_marker_Mm_bm_GSE282765$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_bm_GSE182001.csv")
eos_marker_Mm_il5tg_bm_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "EoP",]
eos_marker_Mm_il5tg_bm_GSE182001 <- eos_marker_Mm_il5tg_bm_GSE182001[order(-eos_marker_Mm_il5tg_bm_GSE182001$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_blood_GSE182001.csv")
eos_marker_Mm_il5tg_blood_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_il5tg_blood_GSE182001 <- eos_marker_Mm_il5tg_blood_GSE182001[order(-eos_marker_Mm_il5tg_blood_GSE182001$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_spleen_GSE182001.csv")
eos_marker_Mm_il5tg_spleen_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_il5tg_spleen_GSE182001 <- eos_marker_Mm_il5tg_spleen_GSE182001[order(-eos_marker_Mm_il5tg_spleen_GSE182001$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_stomach_GSE182001.csv")
eos_marker_Mm_il5tg_stomach_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_il5tg_stomach_GSE182001 <- eos_marker_Mm_il5tg_stomach_GSE182001[order(-eos_marker_Mm_il5tg_stomach_GSE182001$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_duodenum_GSE182001.csv")
eos_marker_Mm_il5tg_duodenum_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_il5tg_duodenum_GSE182001 <- eos_marker_Mm_il5tg_duodenum_GSE182001[order(-eos_marker_Mm_il5tg_duodenum_GSE182001$avg_log2FC), ][1:200, ]
df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_il5tg_colon_GSE182001.csv")
eos_marker_Mm_il5tg_colon_GSE182001 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_il5tg_colon_GSE182001 <- eos_marker_Mm_il5tg_colon_GSE182001[order(-eos_marker_Mm_il5tg_colon_GSE182001$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_bm_GSM7919060.csv")
eos_marker_Mm_bm_GSM7919060 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "EoP",]
eos_marker_Mm_bm_GSM7919060 <- eos_marker_Mm_bm_GSM7919060[order(-eos_marker_Mm_bm_GSM7919060$avg_log2FC), ][1:200, ]

df <- read.csv("/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Mm_Liver_GSE216189.csv")
eos_marker_Mm_Liver_GSE216189 <- df[df$p_val_adj <= 0.05 & df$cluster %in% "Eosinophils",]
eos_marker_Mm_Liver_GSE216189 <- eos_marker_Mm_Liver_GSE216189[order(-eos_marker_Mm_Liver_GSE216189$avg_log2FC), ][1:200, ]

## EoP 
intersectEoP <- Reduce(intersect, list(eos_marker_Mm_bm_GSE282765$gene, eos_marker_Mm_il5tg_bm_GSE182001$gene,eos_marker_Mm_bm_GSM7919060$gene))
x <- list("bm_GSE282765" = eos_marker_Mm_bm_GSE282765$gene, "bm_il5tg_GSE182001" =eos_marker_Mm_il5tg_bm_GSE182001$gene,
          "Mm_bm_GSM7919060" = eos_marker_Mm_bm_GSM7919060$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#76B8E5")

# unique 
setdiff(eos_marker_Mm_bm_GSE282765$gene,Reduce(union, list(eos_marker_Mm_il5tg_bm_GSE182001$gene,eos_marker_Mm_bm_GSM7919060$gene)))
setdiff(eos_marker_Mm_il5tg_bm_GSE182001$gene,Reduce(union, list(eos_marker_Mm_bm_GSE282765$gene,eos_marker_Mm_bm_GSM7919060$gene)))
setdiff(eos_marker_Mm_bm_GSM7919060$gene,Reduce(union, list(eos_marker_Mm_il5tg_bm_GSE182001$gene,eos_marker_Mm_bm_GSE282765$gene)))

# shared 
intersect <- Reduce(intersect, list(eos_marker_Mm_bm_GSE282765$gene,eos_marker_Mm_il5tg_bm_GSE182001$gene))
setdiff(intersect,intersectEoP)

intersect <- Reduce(intersect, list(eos_marker_Mm_bm_GSE282765$gene,eos_marker_Mm_bm_GSM7919060$gene))
setdiff(intersect,intersectEoP)

intersect <- Reduce(intersect, list(eos_marker_Mm_bm_GSM7919060$gene,eos_marker_Mm_il5tg_bm_GSE182001$gene))
setdiff(intersect,intersectEoP)

## Eosinophils lymphatic system 
intersectLymph <- Reduce(intersect, list(eos_marker_Mm_blood_GSE282765$gene,eos_marker_Mm_il5tg_blood_GSE182001$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene ))
x <- list("blood_GSE282765" = eos_marker_Mm_blood_GSE282765$gene, "Mm_il5tg_spleen_GSE182001" =eos_marker_Mm_il5tg_spleen_GSE182001$gene,
          "Mm_il5tg_blood_GSE182001" = eos_marker_Mm_il5tg_blood_GSE182001$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#B5DCB2")

# unique 
setdiff(eos_marker_Mm_blood_GSE282765$gene,Reduce(union, list(eos_marker_Mm_il5tg_blood_GSE182001$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene,intersectLymph)))
setdiff(eos_marker_Mm_il5tg_blood_GSE182001$gene,Reduce(union, list(eos_marker_Mm_blood_GSE282765$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene,intersectLymph)))
setdiff(eos_marker_Mm_il5tg_spleen_GSE182001$gene,Reduce(union, list(eos_marker_Mm_blood_GSE282765$gene,eos_marker_Mm_il5tg_blood_GSE182001$gene),intersectLymph))

# shared 
intersect <- Reduce(intersect, list(eos_marker_Mm_blood_GSE282765$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene))
setdiff(intersect,intersectLymph)

intersect <- Reduce(intersect, list(eos_marker_Mm_il5tg_blood_GSE182001$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene))
setdiff(intersect,intersectLymph)

intersect <- Reduce(intersect, list(eos_marker_Mm_il5tg_blood_GSE182001$gene,eos_marker_Mm_blood_GSE282765$gene))
setdiff(intersect,intersectLymph)

## Eosinophils periphery 
intersect_periphery <- Reduce(intersect, list(eos_marker_Mm_il5tg_stomach_GSE182001$gene, 
                                              eos_marker_Mm_il5tg_duodenum_GSE182001$gene,eos_marker_Mm_il5tg_colon_GSE182001$gene ,
                                              eos_marker_Mm_colon_GSE282765$gene, eos_marker_Mm_tumor_GSE282765$gene, 
                                              eos_marker_Mm_Liver_GSE216189$gene))

intersect_periphery_and_lymphatic <- Reduce(intersect, list(eos_marker_Mm_blood_GSE282765$gene,eos_marker_Mm_il5tg_blood_GSE182001$gene,eos_marker_Mm_il5tg_spleen_GSE182001$gene, 
                                                            eos_marker_Mm_il5tg_stomach_GSE182001$gene, 
                                                            eos_marker_Mm_il5tg_duodenum_GSE182001$gene,eos_marker_Mm_il5tg_colon_GSE182001$gene ,
                                                            eos_marker_Mm_colon_GSE282765$gene, eos_marker_Mm_tumor_GSE282765$gene, 
                                                            eos_marker_Mm_Liver_GSE216189$gene ))

### colon tissues 
x <- list( "colon_GSE282765" = eos_marker_Mm_colon_GSE282765$gene, "Mm_tumor_GSE282765" =eos_marker_Mm_tumor_GSE282765$gene,
           "colon_il5tg_GSE182001" = eos_marker_Mm_il5tg_colon_GSE182001$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#EABBC4")

intersect_colon <- Reduce(intersect, list(eos_marker_Mm_il5tg_colon_GSE182001$gene ,
                                          eos_marker_Mm_colon_GSE282765$gene, eos_marker_Mm_tumor_GSE282765$gene))

# unique 
setdiff(eos_marker_Mm_il5tg_colon_GSE182001$gene,Reduce(union, list(eos_marker_Mm_colon_GSE282765$gene,eos_marker_Mm_tumor_GSE282765$gene)))
setdiff(eos_marker_Mm_colon_GSE282765$gene,Reduce(union, list(eos_marker_Mm_il5tg_colon_GSE182001$gene,eos_marker_Mm_tumor_GSE282765$gene)))
setdiff(eos_marker_Mm_tumor_GSE282765$gene,Reduce(union, list(eos_marker_Mm_il5tg_colon_GSE182001$gene,eos_marker_Mm_colon_GSE282765$gene)))

# shared 
intersect <- Reduce(intersect, list(eos_marker_Mm_il5tg_colon_GSE182001$gene,eos_marker_Mm_colon_GSE282765$gene))
setdiff(intersect,intersect_colon)

intersect <- Reduce(intersect, list(eos_marker_Mm_il5tg_colon_GSE182001$gene,eos_marker_Mm_tumor_GSE282765$gene))
setdiff(intersect,intersect_colon)

intersect <- Reduce(intersect, list(eos_marker_Mm_tumor_GSE282765$gene,eos_marker_Mm_colon_GSE282765$gene))
setdiff(intersect,intersect_colon)

### other peripeheral tissue 
x <- list(  "Liver_GSE216189" = eos_marker_Mm_Liver_GSE216189$gene,"stomach_il5tg_GSE182001" = eos_marker_Mm_il5tg_stomach_GSE182001$gene,
            "duodenum_il5tg_GSE182001" = eos_marker_Mm_il5tg_duodenum_GSE182001$gene)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "purple")

intersect_other_periphery <- Reduce(intersect, list(eos_marker_Mm_Liver_GSE216189$gene ,
                                                    eos_marker_Mm_il5tg_stomach_GSE182001$gene, eos_marker_Mm_il5tg_duodenum_GSE182001$gene))

# unique 
setdiff(eos_marker_Mm_Liver_GSE216189$gene,Reduce(union, list(eos_marker_Mm_il5tg_stomach_GSE182001$gene,eos_marker_Mm_il5tg_duodenum_GSE182001$gene)))
setdiff(eos_marker_Mm_il5tg_stomach_GSE182001$gene,Reduce(union, list(eos_marker_Mm_Liver_GSE216189$gene,eos_marker_Mm_il5tg_duodenum_GSE182001$gene)))
setdiff(eos_marker_Mm_il5tg_duodenum_GSE182001$gene,Reduce(union, list(eos_marker_Mm_Liver_GSE216189$gene,eos_marker_Mm_il5tg_stomach_GSE182001$gene)))

# shared 
intersect <- Reduce(intersect, list(eos_marker_Mm_Liver_GSE216189$gene,eos_marker_Mm_il5tg_stomach_GSE182001$gene))
setdiff(intersect,intersect_other_periphery)

intersect <- Reduce(intersect, list(eos_marker_Mm_il5tg_duodenum_GSE182001$gene,eos_marker_Mm_il5tg_stomach_GSE182001$gene))
setdiff(intersect,intersect_other_periphery)

intersect <- Reduce(intersect, list(eos_marker_Mm_Liver_GSE216189$gene,eos_marker_Mm_il5tg_duodenum_GSE182001$gene))
setdiff(intersect,intersect_other_periphery)
