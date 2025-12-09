######### This code compares SCDC prop.est.mvw scores from different cell types across different doublet detection tools ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load seurat objects 
gene_counts_obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated_gene_count_cutoff.rds")
scDblFinder_obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated_scDblFinder0025.rds")
DoubletFinder_obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated_DoubletFinder.rds")
cell_hashing_obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_multiplets_from_cell_hashing.rds")

# Add method ID 
gene_counts_obj$method <- "gene_counts"
scDblFinder_obj$method <- "scDblFinder"
DoubletFinder_obj$method <- "DoubletFinder"
cell_hashing_obj$method <- "Cell_hashing"

##### Load SCDC results
gene_counts_df <- read.csv("/scratch/khandl/technical/figures/Doublet/upperFeature_cutoff_doublets_deconvolution_result.csv")
scDblFinder_df <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_doublets_deconvolution_result0025.csv")
DoubletFinder_df <- read.csv("/scratch/khandl/technical/figures/Doublet/DoubletFinder_wo_doublets_deconvolution_result2.csv")
cell_hashing_df <- read.csv("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublets_deconvolution_result.csv")

### Extract cells that are present in df 
gene_counts_obj <- subset(gene_counts_obj, cells = gene_counts_df$X)
scDblFinder_obj <- subset(scDblFinder_obj, cells = scDblFinder_df$X)
DoubletFinder_obj <- subset(DoubletFinder_obj, cells = DoubletFinder_df$X)
cell_hashing_obj <- subset(cell_hashing_obj, cells = cell_hashing_df$X)

##### Add SCDC scores to each object 
all(Cells(gene_counts_obj) == gene_counts_df$X)  # should return TRUE
gene_counts_obj <- AddMetaData(gene_counts_obj, metadata = gene_counts_df[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                          "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

all(Cells(scDblFinder_obj) == scDblFinder_df$X)  # should return TRUE
scDblFinder_obj <- AddMetaData(scDblFinder_obj, metadata = scDblFinder_df[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                                              "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

all(Cells(DoubletFinder_obj) == DoubletFinder_df$X)  # should return TRUE
DoubletFinder_obj <- AddMetaData(DoubletFinder_obj, metadata = DoubletFinder_df[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                                              "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

all(Cells(cell_hashing_obj) == cell_hashing_df$X)  # should return TRUE
cell_hashing_obj <- AddMetaData(cell_hashing_obj, metadata = cell_hashing_df[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                                                    "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

##### merge object 
obj <- merge(DoubletFinder_obj, c(scDblFinder_obj,gene_counts_obj,cell_hashing_obj))

##### Plot SCDC scores cross methods 
### Eosinophils
# DoubletFinder
obj <- DoubletFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

# Gene counts
obj <- gene_counts_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

# scDblFinder
obj <- scDblFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df3 <- bind_rows(df_list)

# cell hashing
obj <- cell_hashing_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df4 <- bind_rows(df_list)

df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)

df$method <- factor(df$method, levels = c("Cell_hashing", "gene_counts", "scDblFinder","DoubletFinder"))
p <- ggplot(df, aes(x = method, y =  mean_SCDC_eos, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "Cell_hashing" = "#F9084D", "scDblFinder" = "#F9084D","gene_counts" = "#F9084D","DoubletFinder"="#F9084D"))
ggsave("/scratch/khandl/technical/figures/Doublet/SCDC_score_eos_across_methods.svg", width = 10, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(mean_SCDC_eos ~ method, data = df)
summary(anova)
TukeyHSD(anova)

### Neutrophils
# DoubletFinder
obj <- DoubletFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_neutro = mean(Neutrophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

# Gene counts
obj <- gene_counts_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_neutro = mean(Neutrophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

# scDblFinder
obj <- scDblFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_neutro = mean(Neutrophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df3 <- bind_rows(df_list)

# cell hashing
obj <- cell_hashing_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_neutro = mean(Neutrophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df4 <- bind_rows(df_list)

df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)

df$method <- factor(df$method, levels = c("Cell_hashing", "gene_counts", "scDblFinder","DoubletFinder"))
p <- ggplot(df, aes(x = method, y =  mean_SCDC_neutro, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "Cell_hashing" = "#BD0BF7", "scDblFinder" = "#BD0BF7","gene_counts" = "#BD0BF7","DoubletFinder"="#BD0BF7"))
ggsave("/scratch/khandl/technical/figures/Doublet/SCDC_score_neutro_across_methods.svg", width = 10, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(mean_SCDC_neutro ~ method, data = df)
summary(anova)
TukeyHSD(anova)

### Macrophages
# DoubletFinder
obj <- DoubletFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

# Gene counts
obj <- gene_counts_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

# scDblFinder
obj <- scDblFinder_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df3 <- bind_rows(df_list)

# cell hashing
obj <- cell_hashing_obj
sample <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df4 <- bind_rows(df_list)

df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)

df$method <- factor(df$method, levels = c("Cell_hashing", "gene_counts", "scDblFinder","DoubletFinder"))
p <- ggplot(df, aes(x = method, y =  mean_SCDC_macs, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "Cell_hashing" = "#18C106", "scDblFinder" = "#18C106","gene_counts" = "#18C106","DoubletFinder"="#18C106"))
ggsave("/scratch/khandl/technical/figures/Doublet/SCDC_score_mac_across_methods.svg", width = 10, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(mean_SCDC_macs ~ method, data = df)
summary(anova)
TukeyHSD(anova)

