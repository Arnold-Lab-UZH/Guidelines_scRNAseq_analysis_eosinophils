######### This code compares SCDC prop.est.mvw scores from scDBlFinder between RNA corrected and non-corrected matrices ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load seurat objects 
scDblFinder_obj_corr <- readRDS(file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated_scDblFinder075.rds"))
scDblFinder_obj_no_corr <- readRDS(file.path(seurat_objects_dir, "Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_annotated_scDblFinder075.rds"))

# Add method ID 
scDblFinder_obj_corr$method <- "corr"
scDblFinder_obj_no_corr$method <- "no_corr"

##### Load SCDC results
scDblFinder_df_corr <- read.csv(file.path(doublet_tables_dir, "scDblFinder_wo_doublets_deconvolution_result075_RNA_corr.csv"))
scDblFinder_df_no_corr <- read.csv(file.path(doublet_tables_dir, "scDblFinder_wo_doublets_deconvolution_result075_no_corr.csv"))

### Extract cells that are present in df 
scDblFinder_obj_corr <- subset(scDblFinder_obj_corr, cells = scDblFinder_df_corr$X)
scDblFinder_obj_no_corr <- subset(scDblFinder_obj_no_corr, cells = scDblFinder_df_no_corr$X)

##### Add SCDC scores to each object 
all(Cells(scDblFinder_obj_corr) == scDblFinder_df_corr$X)  # should return TRUE
scDblFinder_obj_corr <- AddMetaData(scDblFinder_obj_corr, metadata = scDblFinder_df_corr[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                                                             "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

all(Cells(scDblFinder_obj_no_corr) == scDblFinder_df_no_corr$X)  # should return TRUE
scDblFinder_obj_no_corr <- AddMetaData(scDblFinder_obj_no_corr, metadata = scDblFinder_df_no_corr[, c("Endothelial", "Epithelial","PCs","Macrophages","DCs","Monocytes",
                                                                                                      "Fibroblasts","T","TAMs","B","Mast","Eosinophils","Neutrophils")])

##### merge object 
obj <- merge(scDblFinder_obj_no_corr, c(scDblFinder_obj_corr))

##### Plot SCDC scores cross methods 
### Eosinophils
# scDblFinder_obj_no_corr
obj <- scDblFinder_obj_no_corr
sample <- (as.data.frame(table(obj$experiment)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

# scDblFinder_obj_corr
obj <- scDblFinder_obj_corr
sample <- (as.data.frame(table(obj$experiment)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_eos = mean(Eosinophils, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = method, y =  mean_SCDC_eos, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "corr" = "#F9084D", "no_corr" = "#F9084D"))
ggsave(file.path(doublet_plots_dir, "SCDC_score_eos_between_corrected_and_non_corrected.svg"), width = 6, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(mean_SCDC_eos ~ method, data = df)
summary(anova)
TukeyHSD(anova)

### Macrophages
# scDblFinder_obj_no_corr
obj <- scDblFinder_obj_no_corr
sample <- (as.data.frame(table(obj$experiment)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

# scDblFinder_obj_corr
obj <- scDblFinder_obj_corr
sample <- (as.data.frame(table(obj$experiment)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "experiment"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(method = .data$method) %>%
    summarise(mean_SCDC_macs = mean(Macrophages, na.rm = TRUE)) %>%
    as.data.frame()
  df$sample <- i
  df_list[[i]] <- df
}
df2 <- bind_rows(df_list)

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = method, y =  mean_SCDC_macs, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "corr" = "#18C106", "no_corr" = "#18C106"))
ggsave(file.path(doublet_plots_dir, "SCDC_score_mac_across_methods_corr_no_corr.svg"), width = 6, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(mean_SCDC_macs ~ method, data = df)
summary(anova)
TukeyHSD(anova)

