########## This code compares gene counts across eosinophil subtypes  ##########
### Datasets used: GSE182001

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load annotated Seruat object 
mouse_data_il5tg <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_eos_annotation.rds"))

##### For each dataset generate a dataframe with the median of gene counts per subtype 
### mouse_data_il5tg
obj <- mouse_data_il5tg
sample <- (as.data.frame(table(obj$sample)))$Var1
df_list <- list()
for (i in sample){
  Idents(obj) <- "sample"
  sub <- subset(obj, idents = i)
  df <- sub@meta.data %>%
    group_by(celltype = .data$annotation) %>%
    summarise(median_nFeature = median(nFeature_RNA, na.rm = TRUE)) %>%
    as.data.frame()
  df$tissue <- i
  df$sample <- i
  df_list[[i]] <- df
}
df <- bind_rows(df_list)

#### per cell type 
## mouse il5tg 
p <- ggplot(df, aes(x = reorder(celltype, median_nFeature, FUN = median), y =  median_nFeature, fill = celltype)) + 
  geom_boxplot(outlier.shape = 16) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme_minimal() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("active" = "#E81818", "circulating" = "#EDE20F","basal" = "#10A069", "immature" = "#E88A1A",
                               "progenitor" = "#26DFED"))
ggsave(file.path(gene_counts_plots_dir, "eos_subtype_il5tg.svg"), width = 8, height = 8, plot = p)

## statistical test --> one way ANOVA 
anova <- aov(median_nFeature ~ celltype, data = df)
summary(anova)
TukeyHSD(anova)

