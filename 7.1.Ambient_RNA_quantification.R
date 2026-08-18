########## This code quantifies ambient RNA content based on IGKC/Igkc expression  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor and Mm colon, NAT, tumor, disseminated 

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load annotated R objects 
human_colon_tumor <- readRDS( file.path(seurat_objects_dir, "Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))
mouse_colon_tumor <- readRDS( file.path(seurat_objects_dir, "Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated.rds"))

##### Plot IGKC/Igkc in FeaturePlots 
Idents(human_colon_tumor) <- "condition"
sub <- subset(human_colon_tumor, idents = c("P1_tumor","P2_tumor","P3_tumor"))
p <- FeaturePlot(sub, features = c("IGKC")) +scale_color_gradientn( colours = c('grey', 'darkred'))
ggsave(file.path(ambient_rna_plots_dir, "IGKC_human.svg"), width = 10, height = 8, plot = p)
# Check where PCs and eosinophils are 
DimPlot(sub, group.by = "annotation", label = TRUE)

Idents(mouse_colon_tumor) <- "condition"
sub <- subset(mouse_colon_tumor, idents = "tumor_wt")
p <- FeaturePlot(sub, features = c("Igkc")) +scale_color_gradientn( colours = c('grey', 'darkred'),limits = c(0,8.5))
ggsave(file.path(ambient_rna_plots_dir, "IGKC_mouse.svg"), width = 10, height = 8, plot = p)
# Check where PCs and eosinophils are 
DimPlot(sub, group.by = "annotation", label = TRUE)

##### Average expression of IGKC/Igkc in Boxplot 
### human_colon 
obj <- human_colon_tumor
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  # Extract IGKC expression + metadata
  dat <- FetchData(sub, vars = c("IGKC", "annotation"), layer = "data")
  
  # Compute average (or median) expression per annotation
  df <- dat %>%
    group_by(annotation) %>%
    summarise(avg_IGKC = median(IGKC, na.rm = TRUE)) %>%
    as.data.frame()
  
  df$sample <- i
  df_list[[i]] <- df
}
df1 <- bind_rows(df_list)

### mouse_colon
obj <- mouse_colon_tumor
conditions <- (as.data.frame(table(obj$condition)))$Var1
df_list <- list()
for (i in conditions){
  Idents(obj) <- "condition"
  sub <- subset(obj, idents = i)
  
  # Extract Igkc expression + metadata
  dat <- FetchData(sub, vars = c("Igkc", "annotation"), layer = "data")
  
  # Compute average (or median) expression per annotation
  df <- dat %>%
    group_by(annotation) %>%
    summarise(avg_IGKC = median(Igkc, na.rm = TRUE)) %>%
    as.data.frame()
  
  df$sample <- i
  df_list[[i]] <- df
  
}
df2 <- bind_rows(df_list)

df_list_all <- list(df1,df2)
df <- bind_rows(df_list_all)

df <- df[df$annotation %in% c("Epithelial","Fibroblasts","Macrophages", "DCs","Endothelial","Monocytes","TAMs","Mast","T","Neutrophils","B","Eosinophils",
                              "PCs"),]

##### Plot boxplot 
p <- ggplot(df, aes(x = reorder(annotation, avg_IGKC, FUN = median), y =  avg_IGKC, fill = annotation)) + 
  geom_boxplot(outlier.shape = 16) + 
  theme_minimal() +   
  #geom_point(position = position_dodge(width = 0.75), size = 2.5,shape = 21, fill = "white",colour = "black", alpha = 0.6, stroke = 1.2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c( "B" = "#F3E972","DCs" = "#E43794","Endothelial" = "#A09167", 
                                "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108",
                                 "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave(file.path(ambient_rna_plots_dir, "IGKC_averag_expression_all.svg"), width = 10, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(avg_IGKC ~ annotation, data = df)
summary(anova)
TukeyHSD(anova)




