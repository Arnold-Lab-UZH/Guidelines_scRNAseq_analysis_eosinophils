######## This code compares the number of genes across cell types between intronic+exonic and exonic reads only mapping  ##########
### Datasets used: GSE282765 and GSE182001

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load annotaed objects 
human_colon <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
human_colon_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_exons_only_annotation.rds")

human_blood <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
human_blood_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_exons_only_annotation.rds")

mouse_colon <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds")
mouse_colon_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_exons_only_annotation.rds")

mouse_blood <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_with_intronic_reads_annotated.rds")
mouse_blood_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_blood_bm_tumor_healthy_forced_cell_determination_exons_only_annotation.rds")

mouse_il5tg <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_with_intronic_reads_annotated.rds")
mouse_il5tg_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_steady_state_forced_cell_determination_exons_only_annotation.rds")

##### Generate a df for plotting with the median number of genes per cell type and smaple 
### human_colon
obj <- human_colon
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Hs"
  df$condition <- i
  df$reads <- "intronic_and_exonic"
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

obj <- human_colon_exon_only
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Hs"
  df$condition <- i
  df$reads <- "exonic_only"
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

### human_blood
obj <- human_blood
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Hs"
  df$condition <- i
  df$reads <- "intronic_and_exonic"
  df_list[[i]] <- df
}
df3 <- bind_rows(df_list)

obj <- human_blood_exon_only
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Hs"
  df$condition <- i
  df$reads <- "exonic_only"
  df_list[[i]] <- df
}
df4 <- bind_rows(df_list)

### mouse_colon
obj <- mouse_colon
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm"
  df$condition <- i
  df$reads <- "intronic_and_exonic"
  df_list[[i]] <- df
}
df5 <- bind_rows(df_list)

obj <- mouse_colon_exon_only
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm"
  df$condition <- i
  df$reads <- "exonic_only"
  df_list[[i]] <- df
}
df6 <- bind_rows(df_list)

### mouse_blood
obj <- mouse_blood
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm"
  df$condition <- i
  df$reads <- "intronic_and_exonic"
  df_list[[i]] <- df
}
df7 <- bind_rows(df_list)

obj <- mouse_blood_exon_only
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm"
  df$condition <- i
  df$reads <- "exonic_only"
  df_list[[i]] <- df
}
df8 <- bind_rows(df_list)

### mouse_il5tg
obj <- mouse_il5tg
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm_il5tg"
  df$condition <- i
  df$reads <- "intronic_and_exonic"
  df_list[[i]] <- df
}
df9 <- bind_rows(df_list)

obj <- mouse_il5tg_exon_only
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$species <- "Mm_il5tg"
  df$condition <- i
  df$reads <- "exonic_only"
  df_list[[i]] <- df
}
df10 <- bind_rows(df_list)

# Combine all df  
df_list <- list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)
df_combined <- bind_rows(df_list)

### Remove ? 
df <- df_combined[df_combined$celltype %in% c("B","Basophils","DCs","Endothelial","EoP","Eosinophils","Epithelial","Fibroblasts","GMPs","HSCs","Macrophages",
                                              "Mast","Monocytes","Neutrophils","PCs","ProMono","ProNeutro","T","TAMs"),]

### Plot in boxplot 
p <- ggplot(df, aes(x=reorder(celltype, median_nFeature, FUN = median), y=as.numeric(median_nFeature), fill=reads)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() +   
  #geom_point(position = position_dodge(width = 0.75), size = 1.5, alpha = 0.6) + 
  scale_fill_manual(values = c("intronic_and_exonic" = "#E87605", "exonic_only" = "#0F92F4")) +
  stat_compare_means(method = "wilcox.test", size = 3) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/intron_exon/all_cell_intron_exon_genes.svg", width = 20, height = 8, plot = p)

### Plot only Eosinophils and EoP
df2 <- df[df$celltype %in% c("Eosinophils","EoP"),]
p <- ggplot(df2, aes(x=reorder(celltype, median_nFeature, FUN = median), y=as.numeric(median_nFeature), fill=reads)) + 
  geom_boxplot(outlier.shape = NA) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  scale_fill_manual(values = c("intronic_and_exonic" = "#E87605", "exonic_only" = "#0F92F4")) +
  stat_compare_means(method = "wilcox.test", size = 3) + theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/intron_exon/Eos_cell_intron_exon_genes.svg", width = 8, height = 8, plot = p)
