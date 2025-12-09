########### This code calculates the gene dropout rate across cell types  ##########
### Datasets used: GSE256088, GSE175930, E-MTAB-14010, GSM7919060, GSE276583, GSE216189, GSE282765, GSE182001

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Technical/1.1.Packages.R")

##### Load R objects 
In_house_data1 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data2 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data3 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data4 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data5 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")

Public_data1 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_GSE256088_et_al_anno.rds")
Public_data2 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_E-MTAB-14010_anno.rds")
Public_data3 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_bm_GSM7919060_anno.rds")
Public_data4 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_blood_GSE276583_anno.rds")
Public_data5 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_liver_GSE216189.rds")
Public_data6 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_esophagus_duodenum_EoE_GSE175930_anno.rds")

##### Remove ? and lowQ and mixed cells 
Idents(In_house_data1) <- "annotation"
In_house_data1 <- subset(In_house_data1, idents = c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                                                    "Monocytes","Neutrophils", "PCs","T","TAMs"))

Idents(In_house_data2) <- "annotation"
In_house_data2 <- subset(In_house_data2, idents = c("Basophils","DCs","Eosinophils","Fibroblasts","Macrophages","Mast",
                                                    "Monocytes","Neutrophils", "PCs","ProNeutro","T"))

Idents(In_house_data3) <- "annotation"
In_house_data3 <- subset(In_house_data3, idents = c("B","DCs","Eosinophils","Macrophages","Mast",
                                                    "Monocytes","Neutrophils", "PCs","T","TAMs"))

Idents(In_house_data4) <- "annotation"
In_house_data4 <- subset(In_house_data4, idents = c("Basophils", "DCs","EoP","Eosinophils", "GMPs","Macrophages","HSCs",
                                                    "Monocytes","Neutrophils","ProMono","ProNeutro","T"))

Idents(In_house_data5) <- "annotation"
In_house_data5 <- subset(In_house_data5, idents = c("Endothelial","EoP", "Eosinophils","Epithelial","Fibroblasts",
                                                    "Macrophages","Neutrophils","PCs","ProNeutro"))

Idents(Public_data1) <- "annotation"
Public_data1 <- subset(Public_data1, idents = c("Neutrophils","T","Monocytes", "Eosinophils","B"))

Idents(Public_data2) <- "annotation"
Public_data2 <- subset(Public_data2, idents = c( "Eosinophils","ProNeutro","Neutrophils"))

Idents(Public_data3) <- "annotation"
Public_data3 <- subset(Public_data3, idents = c("ProNeutro","EoP","GMPs","Neutrophils"))

Idents(Public_data4) <- "annotation"
Public_data4 <- subset(Public_data4, idents = c("Eosinophils","T","Monocytes", "PCs","ProNeutro", "Basophils","B","Mast", "DCs","Neutrophils"))

Idents(Public_data5) <- "annotation"
Public_data5 <- subset(Public_data5, idents = c( "B","Basophils", "DCs","Endothelial","Eosinophils","Epithelial","Hepatocytes",
                                                 "Kupffer","Macrophages","Monocytes","Neutrophils","PCs","Stellate","T"))

Idents(Public_data6) <- "annotation"
Public_data6 <- subset(Public_data6, idents = c("T", "Mast","PCs","Macrophages","Epithelial","Eosinophils","Fibroblasts", "B"))

##### For each dataset generate a dataframe that calculates the dropout rate 
### In_house_data1
obj <- In_house_data1
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "colon")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "colon"
  df$technology <- "BD"
  df$species <- "Hs"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

obj <- In_house_data1
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "tumor")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "tumor"
  df$technology <- "BD"
  df$species <- "Hs"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

### In_house_data2
obj <- In_house_data2
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "blood"
  df$technology <- "BD"
  df$species <- "Hs"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df
}
df3 <- bind_rows(df_list)

### In_house_data3
obj <- In_house_data3
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "colon")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "colon"
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df
}
df4 <- bind_rows(df_list)

obj <- In_house_data3
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "tumor")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "tumor"
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df
}
df5 <- bind_rows(df_list)

### In_house_data4
obj <- In_house_data4
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "blood")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "blood"
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df
}
df6 <- bind_rows(df_list)

obj <- In_house_data4
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "bm")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "bm"
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df
}
df7 <- bind_rows(df_list)

### In_house_data5
obj <- In_house_data5
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- i
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSE282765_Mm_il5tg"
  df_list[[i]] <- df
}
df8 <- bind_rows(df_list)

### Public_data1
obj <- Public_data1
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "blood"
  df$technology <- "HIVE"
  df$species <- "Hs"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSE256088_Hs"
  df_list[[i]] <- df
}
df9 <- bind_rows(df_list)

### Public_data2
obj <- Public_data2
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "blood"
  df$technology <- "BD"
  df$species <- "Hs"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "E_MTAB_14010_Hs"
  df_list[[i]] <- df
}
df10 <- bind_rows(df_list)

### Public_data3
obj <- Public_data3
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "bm"
  df$technology <- "10X"
  df$species <- "Mm"
  df$enrichment <- "Eosinophils"
  df$sample <- i
  df$dataset <- "GSM7919060_Mm"
  df_list[[i]] <- df
}
df11 <- bind_rows(df_list)

### Public_data4
obj <- Public_data4
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "blood"
  df$technology <- "10X"
  df$species <- "Hs"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE276583_Hs"
  df_list[[i]] <- df
}
df12 <- bind_rows(df_list)

### Public_data5
obj <- Public_data5
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "liver"
  df$technology <- "BD"
  df$species <- "Mm"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE216189_Mm"
  df_list[[i]] <- df
}
df13 <- bind_rows(df_list)

### Public_data6
obj <- Public_data6
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "duodenum")
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "duodenum"
  df$technology <- "Seq_Well"
  df$species <- "Hs"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE175930_Hs"
  df_list[[i]] <- df
}
df14 <- bind_rows(df_list)

obj <- Public_data6
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "esophagus")
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  
  df <- data.frame(celltype = unique(sub$annotation)) %>%
    rowwise() %>%
    mutate(
      dropout_rate = (sum(counts[, sub$annotation == celltype, drop = FALSE] == 0) /
                        length(counts[, sub$annotation == celltype, drop = FALSE])) * 100
    ) %>%
    as.data.frame()
  
  df$tissue <- "esophagus"
  df$technology <- "Seq_Well"
  df$species <- "Hs"
  df$enrichment <- "CD45"
  df$sample <- i
  df$dataset <- "GSE175930_Hs"
  df_list[[i]] <- df
}
df15 <- bind_rows(df_list)

#### Combine all dataframes 
df_list_all <- list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15)
df <- bind_rows(df_list_all)

##### Per cell type 
p <- ggplot(df, aes(x = reorder(celltype, dropout_rate, FUN = median), y =  dropout_rate, fill = celltype)) + 
  geom_boxplot(outlier.shape = 16) + 
  theme_minimal() +   
  #geom_point(position = position_dodge(width = 0.75), size = 2.5,shape = 21, fill = "white",colour = "black", alpha = 0.6, stroke = 1.2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c( "B" = "#F3E972","Basophils" = "#EDD6A2", "DCs" = "#E43794","Endothelial" = "#A09167", "EoP" = "#E8AAAA",
                                "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", "GMPs" = "#CC820D",
                                "Hepatocytes"="#EF670A","Kupffer"="#2F5B36", "Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                "HSCs" = "#EFE9BF","Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108",
                                "ProMono"="#D2EFD0", "ProNeutro" = "#D9C1E8", "T" = "#5BC7D9",   "TAMs" = "#516D38", "Stellate"="#877864"))
ggsave("/scratch/khandl/technical/figures/Dropouts/all.svg", width = 20, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(dropout_rate ~ celltype, data = df)
summary(anova)
TukeyHSD(anova)
