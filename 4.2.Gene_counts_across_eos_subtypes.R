########## This code compares the mesians of number of features genes of annotated cell types across datasets  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.0.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Functions.R")

##### read R objects 
mouse_data_il5tg <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_annotation.rds")

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
df_mouse_data_il5tg <- bind_rows(df_list)

#### per cell type 
## mouse il5tg 
p <- ggplot(df_mouse_data_il5tg, aes(x = reorder(celltype, median_nFeature, FUN = median), y =  median_nFeature, fill = celltype)) + 
  geom_boxplot(outlier.shape = 16) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme_minimal() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("active" = "#E81818", "circulating" = "#EDE20F","basal" = "#10A069", "immature" = "#E88A1A",
                               "progenitor" = "#26DFED"))
ggsave("/scratch/khandl/technical/figures/gene_counts_across_datasets/eos_subtype_il5tg.svg", width = 8, height = 8, plot = p)

## statistical test --> one way ANOVA 
anova <- aov(median_nFeature ~ celltype, data = df_mouse_data_il5tg)
summary(anova)
TukeyHSD(anova)

