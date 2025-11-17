########## This code quantifies ambient RNA content using SoupX  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

# GSE282765 
# take forced cell determination and intronic + exonic reads

##### scCDC
### load dataframes 
df1 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_mM.csv")
df2 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_Hs.csv")

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + ylim(0,0.1) + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/scCDC_contamination_per_celltypeMm_Hs.svg", width = 15, height = 8, plot = p)

##### dexontX
### load dataframes 
df1 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/decontX_mM.csv")
df2 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/decontX_Hs.csv")

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(celltype,median_decontX, FUN = median), y =  median_decontX, fill = celltype)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + ylim(0,0.1) + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/decontX_contamination_per_celltypeMm_Hs.svg", width = 15, height = 8, plot = p)

anova <- aov(contamination_values ~ cell_types, data = df)
summary(anova)
TukeyHSD(anova)

##### SoupX
### load dataframes 
df1 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/SoupX_mM.csv")
df2 <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/SoupX_Hs.csv")

df <- rbind(df1,df2)

p <- ggplot(df, aes(x = reorder(cell_types,contamination_values, FUN = median), y =  contamination_values, fill = cell_types)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + scale_y_break(c(35, 75))   + 
  #geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_manual(values =  c( "B" = "#F3E972", "DCs" = "#E43794","Endothelial" = "#A09167", 
                                 "Eosinophils" = "#E22F27", "Epithelial" = "#6D5421",  "Fibroblasts" = "#443511", 
                                 "Hepatocytes"="#EF670A","Macrophages" = "#82C341","Mast" = "#7F7F79", 
                                 "Monocytes" = "#ADD8AB", "Neutrophils" = "#9518ED",  "PCs" = "#B4C108", "T" = "#5BC7D9",   "TAMs" = "#516D38"))
ggsave("/scratch/khandl/technical/figures/Ambient_RNA/SoupX_contamination_per_celltypeMm_Hs.svg", width = 10, height = 8, plot = p)

