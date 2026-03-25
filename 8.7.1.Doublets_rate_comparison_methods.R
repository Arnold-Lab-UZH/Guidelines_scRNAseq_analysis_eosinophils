######### This code compares doublet rates between different doublet detection tools ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load data frames from each algorithm 
### True control 
df1 <- read.csv("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublet_rate.csv")

### no RNA correction 
df2 <- read.csv("/scratch/khandl/technical/figures/Doublet/nFeature_upper_cutoff_doublet_rate_no_corr.csv")
df3 <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_doublet_rate075_no_corr.csv")
df4 <- read.csv("/scratch/khandl/technical/figures/Doublet/DoubletFinder_doublet_rate_no_corr.csv")

### RNA correction 
df5 <- read.csv("/scratch/khandl/technical/figures/Doublet/nFeature_upper_cutoff_doublet_rate_corr.csv")
df6 <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_doublet_rate075_RNA_corr.csv")
df7 <- read.csv("/scratch/khandl/technical/figures/Doublet/DoubletFinder_doublet_rate_corr.csv")

df1$method <- "cell_hashing"
df2$method <- "nFeature_upper_cutoff_no_corr"
df3$method <- "scDblFinder_no_corr"
df4$method <- "DoubletFinder_no_corr"
df5$method <- "nFeature_upper_cutoff_corr"
df6$method <- "scDblFinder_corr"
df7$method <- "DoubletFinder_corr"

## Merge dataframes 
df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)
df <- rbind(df,df5)
df <- rbind(df,df6)
df <- rbind(df,df7)

### Plot in bloxplot 
p <- ggplot(df, aes(x = method, y =  multiplet_rate_per_sample, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
print(p)

##### plot only the corrected ones for the manuscript, corr and no corr are very similar 
### True control 
df1 <- read.csv("/scratch/khandl/technical/figures/Doublet/cell_hashing_doublet_rate.csv")

### RNA correction 
df2 <- read.csv("/scratch/khandl/technical/figures/Doublet/nFeature_upper_cutoff_doublet_rate_corr.csv")
df3 <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_doublet_rate075_RNA_corr.csv")
df4 <- read.csv("/scratch/khandl/technical/figures/Doublet/DoubletFinder_doublet_rate_corr.csv")

df1$method <- "cell_hashing"
df2$method <- "nFeature_upper_cutoff"
df3$method <- "scDblFinder"
df4$method <- "DoubletFinder"

## Merge dataframes 
df <- rbind(df1,df2)
df <- rbind(df,df3)
df <- rbind(df,df4)

### Plot in bloxplot 
df$method <- factor(df$method, levels = c("cell_hashing", "nFeature_upper_cutoff", "scDblFinder","DoubletFinder"))
p <- ggplot(df, aes(x = method, y =  multiplet_rate_per_sample, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values =  c( "cell_hashing" = "#BF4BEF", "scDblFinder" = "#0FD367","nFeature_upper_cutoff" = "#175ABC","DoubletFinder"="#EDE60A"))
ggsave("/scratch/khandl/technical/figures/Doublet/Doublet_rate.svg", width = 8, height = 8, plot = p)

## Statistical test --> one way ANOVA 
anova <- aov(multiplet_rate_per_sample ~ method, data = df)
summary(anova)
TukeyHSD(anova)

