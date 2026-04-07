########### This code calculates the technical gene dropout rate across cell types  ##########
### Datasets used: GSE256088, GSE175930, E-MTAB-14010, GSM7919060, GSE276583, GSE216189, GSE282765, GSE182001

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Technical/1.1.Packages.R")

##### Load R objects 
In_house_data1 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data2 <- readRDS("/scratch/khandl/4.Technical/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data3 <- readRDS("/scratch/khandl/4.Technical/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data4 <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
In_house_data5 <- readRDS("/scratch/khandl/4.Technical/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")

Public_data1 <- readRDS("/scratch/khandl/4.Technical/Hs_blood_GSE256088_et_al_anno.rds")
Public_data2 <- readRDS("/scratch/khandl/4.Technical/Hs_blood_E-MTAB-14010_anno.rds")
Public_data3 <- readRDS("/scratch/khandl/4.Technical/Mm_bm_GSM7919060_anno.rds")
Public_data4 <- readRDS("/scratch/khandl/4.Technical/Hs_blood_GSE276583_anno.rds")
Public_data5 <- readRDS("/scratch/khandl/4.Technical/Mm_liver_GSE216189.rds")
Public_data6 <- readRDS("/scratch/khandl/4.Technical/Hs_esophagus_duodenum_EoE_GSE175930_anno.rds")

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
sample <- c("P1_tissue_ctrl", "P3_tissue_ctrl", "P4_tissue_ctrl", "P5_tissue_ctrl", "P6_tissue_ctrl")
df_list <- list()
for (i in sample){
  # subset sample 
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  # generate list for cell types 
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  # loop through cell types 
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    # extract counts 
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    # skip if cell type not present in sample 
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "colon"
  df1$technology <- "BD"
  df1$species <- "Hs"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df1
}
df2 <- bind_rows(df_list)

obj <- In_house_data1
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "tumor")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "tumor"
  df1$technology <- "BD"
  df1$species <- "Hs"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df1
}
df3 <- bind_rows(df_list)

### In_house_data2
obj <- In_house_data2
sample <- (as.data.frame(table(obj$sample)))$Var1
sample <- c("H3_blood","H5_blood","H7_blood","H8_blood","H9_blood", "P1_blood", "P6_blood", "P7_blood")
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "blood"
  df1$technology <- "BD"
  df1$species <- "Hs"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSE282765_Hs"
  df_list[[i]] <- df1
}
df4 <- bind_rows(df_list)

### In_house_data3
obj <- In_house_data3
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "colon")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "colon"
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df1
}
df5 <- bind_rows(df_list)

obj <- In_house_data3
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "tumor")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "tumor"
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df1
}
df6 <- bind_rows(df_list)

### In_house_data4
obj <- In_house_data4
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "blood")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "blood"
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df1
}
df7 <- bind_rows(df_list)

obj <- In_house_data4
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "bm")
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "bm"
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSE282765_Mm_wt"
  df_list[[i]] <- df1
}
df8 <- bind_rows(df_list)

### In_house_data5
obj <- In_house_data5
sample <- c( "blood", "small_int", "spleen",    "stomach" )
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- i
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSE282765_Mm_il5tg"
  df_list[[i]] <- df1
}
df9 <- bind_rows(df_list)

### Public_data1
obj <- Public_data1
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "blood"
  df1$technology <- "HIVE"
  df1$species <- "Hs"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSE256088_Hs"
  df_list[[i]] <- df1
}
df10 <- bind_rows(df_list)

### Public_data2
obj <- Public_data2
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "blood"
  df1$technology <- "BD"
  df1$species <- "Hs"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "E_MTAB_14010_Hs"
  df_list[[i]] <- df1
}
df11 <- bind_rows(df_list)

### Public_data3
obj <- Public_data3
obj$sample <- obj$condition
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "bm"
  df1$technology <- "10X"
  df1$species <- "Mm"
  df1$enrichment <- "Eosinophils"
  df1$sample <- i
  df1$dataset <- "GSM7919060_Mm"
  df_list[[i]] <- df1
}
df12 <- bind_rows(df_list)

### Public_data4
obj <- Public_data4
obj$sample <- obj$condition
sample <- c("healthy_MNC1", "healthy_MNC2", "healthy_MNC3",  "healthy_PB2", "healthy_PB4", "mild_asthma_PB1","severe_asthma_PB1")
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]

    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "blood"
  df1$technology <- "10X"
  df1$species <- "Hs"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE276583_Hs"
  df_list[[i]] <- df1
}
df13 <- bind_rows(df_list)

### Public_data5
obj <- Public_data5
obj$sample <- obj$condition
sample <- c("Liver_healthy1", "Liver_mets1", "Liver_mets10",   "Liver_mets11",    "Liver_mets3",    "Liver_mets6",    "Liver_mets8",    "Liver_mets9" )
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "liver"
  df1$technology <- "BD"
  df1$species <- "Mm"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE216189_Mm"
  df_list[[i]] <- df1
}
df14 <- bind_rows(df_list)

### Public_data6
obj <- Public_data6
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "duodenum")
sample <- c("Duodenum_EoE_active1", "Duodenum_EoE_active2", "Duodenum_EoE_active3", "Duodenum_EoE_rem1",    "Duodenum_EoE_rem3" )
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "duodenum"
  df1$technology <- "Seq_Well"
  df1$species <- "Hs"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE175930_Hs"
  df_list[[i]] <- df1
}
df15 <- bind_rows(df_list)

obj <- Public_data6
Idents(obj) <- "tissue"
obj <- subset(obj, idents = "esophagus")
sample <- c( "Esophagus_EoE_active2", "Esophagus_EoE_active3", "Esophagus_EoE_rem1",    "Esophagus_EoE_rem4" )
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  df_list2 <- list()
  cell_types_sample <- unique(sub$annotation)
  for (j in cell_types_sample) {
    Idents(sub) <- "annotation"
    sub2 <- subset(sub, idents = j)
    
    counts <- GetAssayData(sub2, assay = "RNA", layer = "counts")
    
    if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) next
    
    # remove genes with overall 0 counts 
    counts <- counts[rowSums(counts) > 0, ]
    
    # calculate dropout rate per cell 
    dropout_rate_cell <- apply(counts, 2, function(x) {
      (sum(x == 0) / length(x)) * 100})
    
    # generate dataframe with cell type and cell ID 
    df_cells <- data.frame(cell = colnames(counts),dropout_rate = dropout_rate_cell)
    
    # calculate the median across cells for each cell type 
    df <- df_cells %>% summarise(median_dropout = median(dropout_rate), .groups = "drop")
    df$cell_type <- j
    df_list2[[j]] <- df
  }
  df1 <- bind_rows(df_list2)
  df1$tissue <- "esophagus"
  df1$technology <- "Seq_Well"
  df1$species <- "Hs"
  df1$enrichment <- "CD45"
  df1$sample <- i
  df1$dataset <- "GSE175930_Hs"
  df_list[[i]] <- df1
}
df16 <- bind_rows(df_list)

#### Combine all dataframes 
df_list_all <- list(df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16)
df <- bind_rows(df_list_all)
length(unique(df$sample)) #69

##### Per cell type 
p <- ggplot(df, aes(x = reorder(cell_type, median_dropout, FUN = median), y =  median_dropout, fill = cell_type)) + 
  geom_boxplot(outlier.shape = 16) + 
  theme_minimal() +   
  #geom_point(position = position_dodge(width = 0.75), size = 2.5,shape = 21, fill = "white",colour = "black", alpha = 0.6, stroke = 1.2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c( "B" = "#F3E972","Basophils" = "#EDD6A2", "DCs" = "#E43794","Endothelial" = "#A09167", "EoP" = "#E8AAAA",
                                "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", "GMPs" = "#CC820D",
                                "Hepatocytes"="#EF670A","Kupffer"="#2F5B36", "Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                "HSCs" = "#EFE9BF","Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108",
                                "ProMono"="#D2EFD0", "ProNeutro" = "#D9C1E8", "T" = "#5BC7D9",   "TAMs" = "#516D38", "Stellate"="#877864"))
ggsave("/scratch/khandl/technical/figures/Dropouts/Dropout_rate_technical.svg", width = 20, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(median_dropout ~ cell_type, data = df)
summary(anova)
TukeyHSD(anova)
