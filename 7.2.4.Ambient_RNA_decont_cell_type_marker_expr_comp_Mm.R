########## This code tests the quality of decontamination tools by their influence on marker gene expression  ##########
### Datasets used: GSE282765; Mm colon healthy, CRC tumor, NAT, disseminated 

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.5.Functions_ambient_RNA_plotting.R")

##### Load R objects 
obj_RNA_scCDC <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_scCDC.rds")
obj_decontX <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_decontX.rds")
obj_SoupX <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_SoupX.rds")

Idents(obj_RNA_scCDC) <- "annotation"
obj_RNA_scCDC <- subset(obj_RNA_scCDC, idents = c("PCs","TAMs","Neutrophils","Eosinophils"))
Idents(obj_decontX) <- "annotation"
obj_decontX <- subset(obj_decontX, idents = c("PCs","TAMs","Neutrophils","Eosinophils"))
Idents(obj_SoupX) <- "annotation"
obj_SoupX <- subset(obj_SoupX, idents = c("PCs","TAMs","Neutrophils","Eosinophils"))

obj_RNA_scCDC <- JoinLayers(obj_RNA_scCDC, assay = "RNA")
obj_RNA_scCDC <- JoinLayers(obj_RNA_scCDC, assay = "Corrected")
obj_decontX <- JoinLayers(obj_decontX)
obj_SoupX <- JoinLayers(obj_SoupX)

### RNA assay 
DefaultAssay(obj_RNA_scCDC) <- "RNA"
sce <- as.SingleCellExperiment(obj_RNA_scCDC)

markers <- list(PC_Markers = c("Igkc","Jchain","Iglc2"),
                TAMs_Marker = c("C1qc","Spp1","Fn1"),
                Neutro_Markers = c("S100a8", "S100a9","Csf3r"),
                Eos_Marker = c("Ccr3","F5","Syne1"))

colData(sce)$cluster <- sce$annotation
df_RNA <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")

### decontX assay 
DefaultAssay(obj_decontX) <- "RNA"
sce <- as.SingleCellExperiment(obj_decontX)

markers <- list(PC_Markers = c("Igkc","Jchain","Iglc2"),
                TAMs_Marker = c("C1qc","Spp1","Fn1"),
                Neutro_Markers = c("S100a8", "S100a9","Csf3r"),
                Eos_Marker = c("Ccr3","F5","Syne1"))

colData(sce)$cluster <- sce$annotation
df_decontX <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")
df_decontX$assay <- "decontX"

### SoupX assay 
DefaultAssay(obj_SoupX) <- "RNA"
sce <- as.SingleCellExperiment(obj_SoupX)

markers <- list(PC_Markers = c("Igkc","Jchain","Iglc2"),
                TAMs_Marker = c("C1qc","Spp1","Fn1"),
                Neutro_Markers = c("S100a8", "S100a9","Csf3r"),
                Eos_Marker = c("Ccr3","F5","Syne1"))

colData(sce)$cluster <- sce$annotation
df_SoupX <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")
df_SoupX$assay <- "SoupX"

### scCDC assay 
DefaultAssay(obj_RNA_scCDC) <- "Corrected"
sce <- as.SingleCellExperiment(obj_RNA_scCDC)

markers <- list(PC_Markers = c("Igkc","Jchain","Iglc2"),
                TAMs_Marker = c("C1qc","Spp1","Fn1"),
                Neutro_Markers = c("S100a8", "S100a9","Csf3r"),
                Eos_Marker = c("Ccr3","F5","Syne1"))

colData(sce)$cluster <- sce$annotation
df_scCDC <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = c("counts"), z = "annotation")
df_scCDC$assay <- "scCDC"

df <- rbind(df_RNA, df_decontX)
df <- rbind(df, df_SoupX)
df <- rbind(df, df_scCDC)

## Plot marker labels of PC markers 
df1 <- df[df$markerLabels %in% "PC_Markers",]

# Set the desired order of assays
df1$assay <- factor(df1$assay, levels = c("counts","decontX","SoupX","scCDC"))

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","decontX"="#0A9B29","scCDC"="#13A2E2","SoupX"="#9C06B7"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/PCmarkersMm.svg", width = 10, height = 8, plot = p)

## Plot marker labels of TAMs markers 
df1 <- df[df$markerLabels %in% "TAMs_Marker",]

# Set the desired order of assays
df1$assay <- factor(df1$assay, levels = c("counts","decontX","SoupX","scCDC"))

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","decontX"="#0A9B29","scCDC"="#13A2E2","SoupX"="#9C06B7"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) + ylim(0,100)
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/TAMsmarkersMm.svg", width = 10, height = 8, plot = p)

## Plot marker labels of Neutrophil markers 
df1 <- df[df$markerLabels %in% "Neutro_Markers",]

# Set the desired order of assays
df1$assay <- factor(df1$assay, levels = c("counts","decontX","SoupX","scCDC"))

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","decontX"="#0A9B29","scCDC"="#13A2E2","SoupX"="#9C06B7"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/NeutromarkersMm.svg", width = 10, height = 8, plot = p)

## Plot marker labels of Eosinophil markers 
df1 <- df[df$markerLabels %in% "Eos_Marker",]

# Set the desired order of assays
df1$assay <- factor(df1$assay, levels = c("counts","decontX","SoupX","scCDC"))

p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#EFD90F","decontX"="#0A9B29","scCDC"="#13A2E2","SoupX"="#9C06B7"), name = "Assay") +  # manual colors & legend title
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +  ylim(0,100)
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/EosmarkersMm.svg", width = 10, height = 8, plot = p)
