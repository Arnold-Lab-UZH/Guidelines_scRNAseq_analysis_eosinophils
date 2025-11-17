########## This code compares outputs of empty droplet characterisation from forced cell determination, EmptyDrops and MALAT1 ##########

# Hs tumor and NAT GSE282765
# take forced cell determination and intronic + exonic reads

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### load objects 
obj_forced <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
obj_automatic <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_automatic_cell_determination_with_intronic_reads_annotated.rds")
obj_emptyDrops <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_emptyDrops_determination_with_intronic_reads_annotated.rds")
obj_Malat1 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_MALAT1_determination_with_intronic_reads_annotated.rds")

## rename ? to Undefined to match the others 
current.cluster.ids <- c("?", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                         "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
new.cluster.ids <- c("Undefined", "B","DCs", "Endothelial","Eosinophils","Epithelial","Fibroblasts",
                     "lowQ", "Macrophages","Mast","Mixed","Monocytes","Neutrophils", "PCs","T","TAMs")
obj_forced$annotation <- plyr::mapvalues(x = obj_forced$annotation, from = current.cluster.ids, to = new.cluster.ids)

##### cell numbers in boxplots 
### Forced BD appraoch
obj <- obj_forced
conditions <- (as.data.frame(table(obj$condition)))$Var1

df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  cond_df <- as.data.frame(table(sub$annotation))
  colnames(cond_df) <- c("annotation",i)
  df_list[[i]] <- cond_df
}
result <- purrr::reduce(df_list, full_join, by = "annotation")
result[is.na(result)] = 0
result <- as.data.frame(result %>%pivot_longer(cols = -annotation,names_to = "condition",values_to = "value"))
results_forced <- result
results_forced$cell_determination <- "BD_forced"

### Automatic BD appraoch
obj <- obj_automatic
conditions <- (as.data.frame(table(obj$condition)))$Var1

df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  cond_df <- as.data.frame(table(sub$annotation))
  colnames(cond_df) <- c("annotation",i)
  df_list[[i]] <- cond_df
}
result <- purrr::reduce(df_list, full_join, by = "annotation")
result[is.na(result)] = 0
result <- as.data.frame(result %>%pivot_longer(cols = -annotation,names_to = "condition",values_to = "value"))
results_automatic <- result
results_automatic$cell_determination <- "BD_automatic"

### obj_emptyDrops appraoch
obj <- obj_emptyDrops
conditions <- (as.data.frame(table(obj$condition)))$Var1

df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  cond_df <- as.data.frame(table(sub$annotation))
  colnames(cond_df) <- c("annotation",i)
  df_list[[i]] <- cond_df
}
result <- purrr::reduce(df_list, full_join, by = "annotation")
result[is.na(result)] = 0
result <- as.data.frame(result %>%pivot_longer(cols = -annotation,names_to = "condition",values_to = "value"))
results_emptyDrops <- result
results_emptyDrops$cell_determination <- "EmptyDrops"

### obj_Malat1 appraoch
obj <- obj_Malat1
conditions <- (as.data.frame(table(obj$condition)))$Var1

df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  cond_df <- as.data.frame(table(sub$annotation))
  colnames(cond_df) <- c("annotation",i)
  df_list[[i]] <- cond_df
}
result <- purrr::reduce(df_list, full_join, by = "annotation")
result[is.na(result)] = 0
result <- as.data.frame(result %>%pivot_longer(cols = -annotation,names_to = "condition",values_to = "value"))
results_Malat1 <- result
results_Malat1$cell_determination <- "MALAT1"

#### combine all dataframes 
df <- rbind(results_forced,results_automatic)
df <- rbind(df, results_emptyDrops)
df <- rbind(df, results_Malat1)

##### remove ? 
df <- df[df$annotation %in% c("B","DCs","Endothelial","Eosinophils","Epithelial","Fibroblasts","lowQ", "Macrophages",
                              "Mast","Monocytes","Neutrophils","PCs","T","TAMs"),]

##### plot in boxplot 
p <- ggplot(df, aes(x=reorder(annotation, value, FUN = median), y=as.numeric(value), fill=cell_determination)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  scale_fill_manual(values = c("BD_forced" = "#F4C80F", "BD_automatic" = "#890FF4","EmptyDrops" ="#60A4F4", "MALAT1"="#1C7F08" ))  +
  theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/summary_empty/across_cell_types.svg", width = 20, height = 8, plot = p)

### pairwise wilcoxon test automatic vs. the rest 
df2 <- df[df$cell_determination %in% c("BD_automatic","BD_forced"),]
ggplot(df2, aes(x=reorder(annotation, value, FUN = median), y=as.numeric(value), fill=cell_determination)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  scale_fill_manual(values = c("BD_forced" = "#F4C80F", "BD_automatic" = "#890FF4","EmptyDrops" ="#60A4F4", "MALAT1"="#1C7F08" ))  +
  theme(axis.text.x = element_text(angle = 45)) +
  stat_compare_means(method = "wilcox.test", size = 3) + theme(axis.text.x = element_text(angle = 45)) 

##### EoP ad Eosinophils in specific 
df2 <- df[df$annotation %in% c("Eosinophils"),]
p <- ggplot(df2, aes(x=reorder(cell_determination, value, FUN = median), y=as.numeric(value), fill=cell_determination)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  scale_fill_manual(values = c("BD_forced" = "#F4C80F", "BD_automatic" = "#890FF4","EmptyDrops" ="#60A4F4", "MALAT1"="#1C7F08" )) + 
   theme(axis.text.x = element_text(angle = 45)) 
ggsave("/scratch/khandl/technical/figures/summary_empty/Eos.svg", width = 8, height = 8, plot = p)

