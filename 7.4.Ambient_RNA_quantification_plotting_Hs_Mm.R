########## This plots contamination levels from different Ambient RNA detection tools ##########
### Datasets used: GSE282765; Mm colon healthy, CRC tumor, NAT, disseminated;   Hs CRC NAT and tumor 

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### scCDC
df1 <- read.csv(file.path(ambient_rna_tables_dir, "scCDC_mM.csv"))
df2 <- read.csv(file.path(ambient_rna_tables_dir, "scCDC_Hs.csv"))

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + ylim(0,0.1) + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave(file.path(ambient_rna_plots_dir, "scCDC_contamination_per_celltypeMm_Hs.svg"), width = 15, height = 8, plot = p)

anova <- aov(contamination_values ~ cell_types, data = df)
summary(anova)
TukeyHSD(anova)

##### dexontX
df1 <- read.csv(file.path(ambient_rna_tables_dir, "decontX_mM.csv"))
df2 <- read.csv(file.path(ambient_rna_tables_dir, "decontX_Hs.csv"))

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(celltype,median_decontX, FUN = median), y =  median_decontX, fill = celltype)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + ylim(0,0.1) + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave(file.path(ambient_rna_plots_dir, "decontX_contamination_per_celltypeMm_Hs.svg"), width = 15, height = 8, plot = p)

anova <- aov(median_decontX ~ celltype, data = df)
summary(anova)
TukeyHSD(anova)

##### SoupX
df1 <- read.csv(file.path(ambient_rna_tables_dir, "SoupX_mM.csv"))
df2 <- read.csv(file.path(ambient_rna_tables_dir, "SoupX_Hs.csv"))

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + scale_y_break(c(35, 75))   + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave(file.path(ambient_rna_plots_dir, "SoupX_contamination_per_celltypeMm_Hs.svg"), width = 10, height = 8, plot = p)

anova <- aov(contamination_values ~ cell_types, data = df)
summary(anova)
TukeyHSD(anova)
