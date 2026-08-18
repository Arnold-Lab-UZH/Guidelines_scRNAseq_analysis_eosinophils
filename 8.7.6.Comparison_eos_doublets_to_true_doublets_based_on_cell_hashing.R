######### This code analyzes eosinophil contributing doublet identities from doublet identification tools to true doublets from cell hashing ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Eosinophil homo- and heterotypic doublets 
### Extract homo and heterotypic barcodes per experiment 
## Cell Hashing = true doublets 
df <- read.csv(file.path(doublet_tables_dir, "cell_hashing_doublets_deconvolution_result.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_hash <- (df_h1[df_h1$Eosinophils > 0.25 & df_h1$Eosinophils <= 0.75,])$X
h2_hetero_hash <- (df_h2[df_h2$Eosinophils > 0.25 & df_h2$Eosinophils <= 0.75,])$X
h3_hetero_hash <- (df_h3[df_h3$Eosinophils > 0.25 & df_h3$Eosinophils <= 0.75,])$X
h4_hetero_hash <- (df_h4[df_h4$Eosinophils > 0.25 & df_h4$Eosinophils <= 0.75,])$X
h5_hetero_hash <- (df_h5[df_h5$Eosinophils > 0.25 & df_h5$Eosinophils <= 0.75,])$X
h6_hetero_hash <- (df_h6[df_h6$Eosinophils > 0.25 & df_h6$Eosinophils <= 0.75,])$X
h7_hetero_hash <- (df_h7[df_h7$Eosinophils > 0.25 & df_h7$Eosinophils <= 0.75,])$X

h1_homo_hash <- (df_h1[df_h1$Eosinophils > 0.75,])$X
h2_homo_hash <- (df_h2[df_h2$Eosinophils > 0.75,])$X
h3_homo_hash <- (df_h3[df_h3$Eosinophils > 0.75,])$X
h4_homo_hash <- (df_h4[df_h4$Eosinophils > 0.75,])$X
h5_homo_hash <- (df_h5[df_h5$Eosinophils > 0.75,])$X
h6_homo_hash <- (df_h6[df_h6$Eosinophils > 0.75,])$X
h7_homo_hash <- (df_h7[df_h7$Eosinophils > 0.75,])$X

## Based on gene counts 
df <- read.csv(file.path(doublet_tables_dir, "upperFeature_cutoff_doublets_deconvolution_result_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_counts <- (df_h1[df_h1$Eosinophils > 0.25 & df_h1$Eosinophils <= 0.75,])$X
h2_hetero_counts <- (df_h2[df_h2$Eosinophils > 0.25 & df_h2$Eosinophils <= 0.75,])$X
h3_hetero_counts <- (df_h3[df_h3$Eosinophils > 0.25 & df_h3$Eosinophils <= 0.75,])$X
h4_hetero_counts <- (df_h4[df_h4$Eosinophils > 0.25 & df_h4$Eosinophils <= 0.75,])$X
h5_hetero_counts <- (df_h5[df_h5$Eosinophils > 0.25 & df_h5$Eosinophils <= 0.75,])$X
h6_hetero_counts <- (df_h6[df_h6$Eosinophils > 0.25 & df_h6$Eosinophils <= 0.75,])$X
h7_hetero_counts <- (df_h7[df_h7$Eosinophils > 0.25 & df_h7$Eosinophils <= 0.75,])$X

h1_homo_counts <- (df_h1[df_h1$Eosinophils > 0.75,])$X
h2_homo_counts <- (df_h2[df_h2$Eosinophils > 0.75,])$X
h3_homo_counts <- (df_h3[df_h3$Eosinophils > 0.75,])$X
h4_homo_counts <- (df_h4[df_h4$Eosinophils > 0.75,])$X
h5_homo_counts <- (df_h5[df_h5$Eosinophils > 0.75,])$X
h6_homo_counts <- (df_h6[df_h6$Eosinophils > 0.75,])$X
h7_homo_counts <- (df_h7[df_h7$Eosinophils > 0.75,])$X

## Based on scDblFinder
df <- read.csv(file.path(doublet_tables_dir, "scDblFinder_wo_doublets_deconvolution_result075_RNA_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_scDbl <- (df_h1[df_h1$Eosinophils > 0.25 & df_h1$Eosinophils <= 0.75,])$X
h2_hetero_scDbl <- (df_h2[df_h2$Eosinophils > 0.25 & df_h2$Eosinophils <= 0.75,])$X
h3_hetero_scDbl <- (df_h3[df_h3$Eosinophils > 0.25 & df_h3$Eosinophils <= 0.75,])$X
h4_hetero_scDbl <- (df_h4[df_h4$Eosinophils > 0.25 & df_h4$Eosinophils <= 0.75,])$X
h5_hetero_scDbl <- (df_h5[df_h5$Eosinophils > 0.25 & df_h5$Eosinophils <= 0.75,])$X
h6_hetero_scDbl <- (df_h6[df_h6$Eosinophils > 0.25 & df_h6$Eosinophils <= 0.75,])$X
h7_hetero_scDbl <- (df_h7[df_h7$Eosinophils > 0.25 & df_h7$Eosinophils <= 0.75,])$X

h1_homo_scDbl <- (df_h1[df_h1$Eosinophils > 0.75,])$X
h2_homo_scDbl <- (df_h2[df_h2$Eosinophils > 0.75,])$X
h3_homo_scDbl <- (df_h3[df_h3$Eosinophils > 0.75,])$X
h4_homo_scDbl <- (df_h4[df_h4$Eosinophils > 0.75,])$X
h5_homo_scDbl <- (df_h5[df_h5$Eosinophils > 0.75,])$X
h6_homo_scDbl <- (df_h6[df_h6$Eosinophils > 0.75,])$X
h7_homo_scDbl <- (df_h7[df_h7$Eosinophils > 0.75,])$X

## Based on DoubletFinder
df <- read.csv(file.path(doublet_tables_dir, "DoubletFinder_wo_doublets_deconvolution_result_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_DoubletF <- (df_h1[df_h1$Eosinophils > 0.25 & df_h1$Eosinophils <= 0.75,])$X
h2_hetero_DoubletF <- (df_h2[df_h2$Eosinophils > 0.25 & df_h2$Eosinophils <= 0.75,])$X
h3_hetero_DoubletF <- (df_h3[df_h3$Eosinophils > 0.25 & df_h3$Eosinophils <= 0.75,])$X
h4_hetero_DoubletF <- (df_h4[df_h4$Eosinophils > 0.25 & df_h4$Eosinophils <= 0.75,])$X
h5_hetero_DoubletF <- (df_h5[df_h5$Eosinophils > 0.25 & df_h5$Eosinophils <= 0.75,])$X
h6_hetero_DoubletF <- (df_h6[df_h6$Eosinophils > 0.25 & df_h6$Eosinophils <= 0.75,])$X
h7_hetero_DoubletF <- (df_h7[df_h7$Eosinophils > 0.25 & df_h7$Eosinophils <= 0.75,])$X

h1_homo_DoubletF <- (df_h1[df_h1$Eosinophils > 0.75,])$X
h2_homo_DoubletF <- (df_h2[df_h2$Eosinophils > 0.75,])$X
h3_homo_DoubletF <- (df_h3[df_h3$Eosinophils > 0.75,])$X
h4_homo_DoubletF <- (df_h4[df_h4$Eosinophils > 0.75,])$X
h5_homo_DoubletF <- (df_h5[df_h5$Eosinophils > 0.75,])$X
h6_homo_DoubletF <- (df_h6[df_h6$Eosinophils > 0.75,])$X
h7_homo_DoubletF <- (df_h7[df_h7$Eosinophils > 0.75,])$X

## Calculate percentage of homo and hetero eosinophil-contributing doublets to ground truth = cell hashing
## gene counts
# heterotypic 
matching_bc <- length(intersect(h1_hetero_counts,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_counts,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_counts,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_counts,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_counts,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_counts,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_counts,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_counts <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_counts) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_counts$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_counts$identity <- "heterotypic"
percentage_hetero_counts$method <- "gene_counts"

# homotypic 
matching_bc <- length(intersect(h1_homo_counts,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_counts,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_counts,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_counts,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_counts,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_counts,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_counts,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_counts <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_counts) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_counts$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_counts$identity <- "homotypic"
percentage_homo_counts$method <- "gene_counts"

## scDblFinder
# heterotypic 
matching_bc <- length(intersect(h1_hetero_scDbl,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_scDbl,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_scDbl,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_scDbl,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_scDbl,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_scDbl,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_scDbl,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_scDbl <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_scDbl) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_scDbl$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_scDbl$identity <- "heterotypic"
percentage_hetero_scDbl$method <- "scDblFinder"

# homotypic 
matching_bc <- length(intersect(h1_homo_scDbl,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_scDbl,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_scDbl,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_scDbl,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_scDbl,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_scDbl,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_scDbl,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_scDbl <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_scDbl) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_scDbl$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_scDbl$identity <- "homotypic"
percentage_homo_scDbl$method <- "scDblFinder"

## gene counts
# heterotypic 
matching_bc <- length(intersect(h1_hetero_DoubletF,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_DoubletF,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_DoubletF,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_DoubletF,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_DoubletF,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_DoubletF,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_DoubletF,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_DoubletF <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_DoubletF) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_DoubletF$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_DoubletF$identity <- "heterotypic"
percentage_hetero_DoubletF$method <- "DoubletFinder"

# homotypic 
matching_bc <- length(intersect(h1_homo_DoubletF,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_DoubletF,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_DoubletF,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_DoubletF,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_DoubletF,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_DoubletF,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_DoubletF,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_DoubletF <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_DoubletF) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_DoubletF$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_DoubletF$identity <- "homotypic"
percentage_homo_DoubletF$method <- "DoubletFinder"

### combine dataframes 
df <- rbind(percentage_hetero_counts, percentage_hetero_DoubletF)
df <- rbind(df, percentage_hetero_scDbl)
df <- rbind(df, percentage_homo_counts)
df <- rbind(df, percentage_homo_DoubletF)
df <- rbind(df, percentage_homo_scDbl)

df$method <- factor(df$method, levels = c( "gene_counts", "scDblFinder","DoubletFinder"))

# heterotype 
df_het <- df[df$identity %in% "heterotypic",]
p <- ggplot(df_het, aes(x = method, y =  Value, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(1,100) + 
  scale_fill_manual(values =  c( "scDblFinder" = "#F9084D","gene_counts" = "#F9084D","DoubletFinder"="#F9084D"))
ggsave(file.path(doublet_plots_dir, "Percentage_of_true_doublets_hashing_eos_heterotypic.svg"), width = 8, height = 8, plot = p)

# homotypic 
df_homo <- df[df$identity %in% "homotypic",]
p <- ggplot(df_homo, aes(x = method, y =  Value, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(1,100) + 
  scale_fill_manual(values =  c( "scDblFinder" = "#F9084D","gene_counts" = "#F9084D","DoubletFinder"="#F9084D"))
ggsave(file.path(doublet_plots_dir, "Percentage_of_true_doublets_hashing_eos_homotypic.svg"), width = 8, height = 8, plot = p)

## Run statistical analysis for each group 
# Heterotypic 
a <- df[df$identity %in% "heterotypic",]
## statistical test --> one way ANOVA 
anova <- aov(Value ~ method, data = a)
summary(anova)
TukeyHSD(anova)

# homotypic 
a <- df[df$identity %in% "homotypic",]
## statistical test --> one way ANOVA 
anova <- aov(Value ~ method, data = a)
summary(anova)
TukeyHSD(anova)

##### Macrophages homo- and heterotypic doublets 
### Extract homo and heterotypic barcodes per experiment 
## Cell Hashing = true doublets 
df <- read.csv(file.path(doublet_tables_dir, "cell_hashing_doublets_deconvolution_result.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_hash <- (df_h1[df_h1$Macrophages > 0.25 & df_h1$Macrophages <= 0.75,])$X
h2_hetero_hash <- (df_h2[df_h2$Macrophages > 0.25 & df_h2$Macrophages <= 0.75,])$X
h3_hetero_hash <- (df_h3[df_h3$Macrophages > 0.25 & df_h3$Macrophages <= 0.75,])$X
h4_hetero_hash <- (df_h4[df_h4$Macrophages > 0.25 & df_h4$Macrophages <= 0.75,])$X
h5_hetero_hash <- (df_h5[df_h5$Macrophages > 0.25 & df_h5$Macrophages <= 0.75,])$X
h6_hetero_hash <- (df_h6[df_h6$Macrophages > 0.25 & df_h6$Macrophages <= 0.75,])$X
h7_hetero_hash <- (df_h7[df_h7$Macrophages > 0.25 & df_h7$Macrophages <= 0.75,])$X

h1_homo_hash <- (df_h1[df_h1$Macrophages > 0.75,])$X
h2_homo_hash <- (df_h2[df_h2$Macrophages > 0.75,])$X
h3_homo_hash <- (df_h3[df_h3$Macrophages > 0.75,])$X
h4_homo_hash <- (df_h4[df_h4$Macrophages > 0.75,])$X
h5_homo_hash <- (df_h5[df_h5$Macrophages > 0.75,])$X
h6_homo_hash <- (df_h6[df_h6$Macrophages > 0.75,])$X
h7_homo_hash <- (df_h7[df_h7$Macrophages > 0.75,])$X

## Based on gene counts 
df <- read.csv(file.path(doublet_tables_dir, "upperFeature_cutoff_doublets_deconvolution_result_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_counts <- (df_h1[df_h1$Macrophages > 0.25 & df_h1$Macrophages <= 0.75,])$X
h2_hetero_counts <- (df_h2[df_h2$Macrophages > 0.25 & df_h2$Macrophages <= 0.75,])$X
h3_hetero_counts <- (df_h3[df_h3$Macrophages > 0.25 & df_h3$Macrophages <= 0.75,])$X
h4_hetero_counts <- (df_h4[df_h4$Macrophages > 0.25 & df_h4$Macrophages <= 0.75,])$X
h5_hetero_counts <- (df_h5[df_h5$Macrophages > 0.25 & df_h5$Macrophages <= 0.75,])$X
h6_hetero_counts <- (df_h6[df_h6$Macrophages > 0.25 & df_h6$Macrophages <= 0.75,])$X
h7_hetero_counts <- (df_h7[df_h7$Macrophages > 0.25 & df_h7$Macrophages <= 0.75,])$X

h1_homo_counts <- (df_h1[df_h1$Macrophages > 0.75,])$X
h2_homo_counts <- (df_h2[df_h2$Macrophages > 0.75,])$X
h3_homo_counts <- (df_h3[df_h3$Macrophages > 0.75,])$X
h4_homo_counts <- (df_h4[df_h4$Macrophages > 0.75,])$X
h5_homo_counts <- (df_h5[df_h5$Macrophages > 0.75,])$X
h6_homo_counts <- (df_h6[df_h6$Macrophages > 0.75,])$X
h7_homo_counts <- (df_h7[df_h7$Macrophages > 0.75,])$X

## Based on scDblFinder
df <- read.csv(file.path(doublet_tables_dir, "scDblFinder_wo_doublets_deconvolution_result075_RNA_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_scDbl <- (df_h1[df_h1$Macrophages > 0.25 & df_h1$Macrophages <= 0.75,])$X
h2_hetero_scDbl <- (df_h2[df_h2$Macrophages > 0.25 & df_h2$Macrophages <= 0.75,])$X
h3_hetero_scDbl <- (df_h3[df_h3$Macrophages > 0.25 & df_h3$Macrophages <= 0.75,])$X
h4_hetero_scDbl <- (df_h4[df_h4$Macrophages > 0.25 & df_h4$Macrophages <= 0.75,])$X
h5_hetero_scDbl <- (df_h5[df_h5$Macrophages > 0.25 & df_h5$Macrophages <= 0.75,])$X
h6_hetero_scDbl <- (df_h6[df_h6$Macrophages > 0.25 & df_h6$Macrophages <= 0.75,])$X
h7_hetero_scDbl <- (df_h7[df_h7$Macrophages > 0.25 & df_h7$Macrophages <= 0.75,])$X

h1_homo_scDbl <- (df_h1[df_h1$Macrophages > 0.75,])$X
h2_homo_scDbl <- (df_h2[df_h2$Macrophages > 0.75,])$X
h3_homo_scDbl <- (df_h3[df_h3$Macrophages > 0.75,])$X
h4_homo_scDbl <- (df_h4[df_h4$Macrophages > 0.75,])$X
h5_homo_scDbl <- (df_h5[df_h5$Macrophages > 0.75,])$X
h6_homo_scDbl <- (df_h6[df_h6$Macrophages > 0.75,])$X
h7_homo_scDbl <- (df_h7[df_h7$Macrophages > 0.75,])$X

## Based on DoubletFinder
df <- read.csv(file.path(doublet_tables_dir, "DoubletFinder_wo_doublets_deconvolution_result_corr.csv"))
df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

h1_hetero_DoubletF <- (df_h1[df_h1$Macrophages > 0.25 & df_h1$Macrophages <= 0.75,])$X
h2_hetero_DoubletF <- (df_h2[df_h2$Macrophages > 0.25 & df_h2$Macrophages <= 0.75,])$X
h3_hetero_DoubletF <- (df_h3[df_h3$Macrophages > 0.25 & df_h3$Macrophages <= 0.75,])$X
h4_hetero_DoubletF <- (df_h4[df_h4$Macrophages > 0.25 & df_h4$Macrophages <= 0.75,])$X
h5_hetero_DoubletF <- (df_h5[df_h5$Macrophages > 0.25 & df_h5$Macrophages <= 0.75,])$X
h6_hetero_DoubletF <- (df_h6[df_h6$Macrophages > 0.25 & df_h6$Macrophages <= 0.75,])$X
h7_hetero_DoubletF <- (df_h7[df_h7$Macrophages > 0.25 & df_h7$Macrophages <= 0.75,])$X

h1_homo_DoubletF <- (df_h1[df_h1$Macrophages > 0.75,])$X
h2_homo_DoubletF <- (df_h2[df_h2$Macrophages > 0.75,])$X
h3_homo_DoubletF <- (df_h3[df_h3$Macrophages > 0.75,])$X
h4_homo_DoubletF <- (df_h4[df_h4$Macrophages > 0.75,])$X
h5_homo_DoubletF <- (df_h5[df_h5$Macrophages > 0.75,])$X
h6_homo_DoubletF <- (df_h6[df_h6$Macrophages > 0.75,])$X
h7_homo_DoubletF <- (df_h7[df_h7$Macrophages > 0.75,])$X

## Calculate percentage of homo and hetero eosinophil-contributing doublets to ground truth = cell hashing
## gene counts
# heterotypic 
matching_bc <- length(intersect(h1_hetero_counts,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_counts,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_counts,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_counts,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_counts,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_counts,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_counts,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_counts <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_counts) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_counts$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_counts$identity <- "heterotypic"
percentage_hetero_counts$method <- "gene_counts"

# homotypic 
matching_bc <- length(intersect(h1_homo_counts,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_counts,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_counts,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_counts,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_counts,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_counts,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_counts,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_counts <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_counts) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_counts$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_counts$identity <- "homotypic"
percentage_homo_counts$method <- "gene_counts"

## scDblFinder
# heterotypic 
matching_bc <- length(intersect(h1_hetero_scDbl,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_scDbl,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_scDbl,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_scDbl,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_scDbl,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_scDbl,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_scDbl,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_scDbl <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_scDbl) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_scDbl$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_scDbl$identity <- "heterotypic"
percentage_hetero_scDbl$method <- "scDblFinder"

# homotypic 
matching_bc <- length(intersect(h1_homo_scDbl,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_scDbl,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_scDbl,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_scDbl,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_scDbl,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_scDbl,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_scDbl,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_scDbl <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_scDbl) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_scDbl$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_scDbl$identity <- "homotypic"
percentage_homo_scDbl$method <- "scDblFinder"

## gene counts
# heterotypic 
matching_bc <- length(intersect(h1_hetero_DoubletF,h1_hetero_hash))
df1 <- matching_bc *100/ length(h1_hetero_hash)

matching_bc <- length(intersect(h2_hetero_DoubletF,h2_hetero_hash))
df2 <- matching_bc *100/ length(h2_hetero_hash)

matching_bc <- length(intersect(h3_hetero_DoubletF,h3_hetero_hash))
df3 <- matching_bc *100/ length(h3_hetero_hash)

matching_bc <- length(intersect(h4_hetero_DoubletF,h4_hetero_hash))
df4 <- matching_bc *100/ length(h4_hetero_hash)

matching_bc <- length(intersect(h5_hetero_DoubletF,h5_hetero_hash))
df5 <- matching_bc *100/ length(h5_hetero_hash)

matching_bc <- length(intersect(h6_hetero_DoubletF,h6_hetero_hash))
df6 <- matching_bc *100/ length(h6_hetero_hash)

matching_bc <- length(intersect(h7_hetero_DoubletF,h7_hetero_hash))
df7 <- matching_bc *100/ length(h7_hetero_hash)

percentage_hetero_DoubletF <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_hetero_DoubletF) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_hetero_DoubletF$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_hetero_DoubletF$identity <- "heterotypic"
percentage_hetero_DoubletF$method <- "DoubletFinder"

# homotypic 
matching_bc <- length(intersect(h1_homo_DoubletF,h1_homo_hash))
df1 <- matching_bc *100/ length(h1_homo_hash)

matching_bc <- length(intersect(h2_homo_DoubletF,h2_homo_hash))
df2 <- matching_bc *100/ length(h2_homo_hash)

matching_bc <- length(intersect(h3_homo_DoubletF,h3_homo_hash))
df3 <- matching_bc *100/ length(h3_homo_hash)

matching_bc <- length(intersect(h4_homo_DoubletF,h4_homo_hash))
df4 <- matching_bc *100/ length(h4_homo_hash)

matching_bc <- length(intersect(h5_homo_DoubletF,h5_homo_hash))
df5 <- matching_bc *100/ length(h5_homo_hash)

matching_bc <- length(intersect(h6_homo_DoubletF,h6_homo_hash))
df6 <- matching_bc *100/ length(h6_homo_hash)

matching_bc <- length(intersect(h7_homo_DoubletF,h7_homo_hash))
df7 <- matching_bc *100/ length(h7_homo_hash)

percentage_homo_DoubletF <- as.data.frame(c(df1,df2,df3,df4,df5,df6,df7))
colnames(percentage_homo_DoubletF) <- "Value"
# experiment 7 = h6/P6, experiment 8 = H7/P7
percentage_homo_DoubletF$exp <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
percentage_homo_DoubletF$identity <- "homotypic"
percentage_homo_DoubletF$method <- "DoubletFinder"

### combine dataframes 
df <- rbind(percentage_hetero_counts, percentage_hetero_DoubletF)
df <- rbind(df, percentage_hetero_scDbl)
df <- rbind(df, percentage_homo_counts)
df <- rbind(df, percentage_homo_DoubletF)
df <- rbind(df, percentage_homo_scDbl)

df$method <- factor(df$method, levels = c( "gene_counts", "scDblFinder","DoubletFinder"))

# heterotype 
df_het <- df[df$identity %in% "heterotypic",]
p <- ggplot(df_het, aes(x = method, y =  Value, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(1,100) + 
  scale_fill_manual(values =  c( "scDblFinder" = "#37B44A","gene_counts" = "#37B44A","DoubletFinder"="#37B44A"))
ggsave(file.path(doublet_plots_dir, "Percentage_of_true_doublets_hashing_mac_heterotypic.svg"), width = 8, height = 8, plot = p)

# homotypic 
df_homo <- df[df$identity %in% "homotypic",]
p <- ggplot(df_homo, aes(x = method, y =  Value, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(1,100) + 
  scale_fill_manual(values =  c( "scDblFinder" = "#37B44A","gene_counts" = "#37B44A","DoubletFinder"="#37B44A"))
ggsave(file.path(doublet_plots_dir, "Percentage_of_true_doublets_hashing_mac_homotypic.svg"), width = 8, height = 8, plot = p)

## Run statistical analysis for each group 
# Heterotypic 
a <- df[df$identity %in% "heterotypic",]
## statistical test --> one way ANOVA 
anova <- aov(Value ~ method, data = a)
summary(anova)
TukeyHSD(anova)

# homotypic 
a <- df[df$identity %in% "homotypic",]
## statistical test --> one way ANOVA 
anova <- aov(Value ~ method, data = a)
summary(anova)
TukeyHSD(anova)
