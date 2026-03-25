######### This code analyzes eosinophil contributing doublets from scDBlFinder outputs comparing datasets with and without RNA correction (scCDC) ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Compare the percentage of heterotypeic eos, homotypeic eos and other based on prop.est.mvw Eos scores 
### Based on scDblFinder
df <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_doublets_deconvolution_result075_RNA_corr.csv")

df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

all_h1 <- length(rownames(df_h1))
doublets_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.25 & df_h1$Eosinophils <= 0.75 ,]))
homo_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.75  ,]))
other_h1 <- length(rownames(df_h1[df_h1$Eosinophils ==0 &df_h1$Eosinophils <= 0.25  ,]))
doublets_h1 <- (doublets_h1*100)/all_h1
homo_h1 <- (homo_h1*100)/all_h1
other_h1 <- (other_h1*100)/all_h1
all_h1 <- (all_h1*100)/all_h1
h1_df <- data.frame(c("all","doublets","homo","other"), c(all_h1,doublets_h1,homo_h1,other_h1))
colnames(h1_df) <- c("identity","values")
h1_df$exp <- "Exp1"

all_h2 <- length(rownames(df_h2))
doublets_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.25 & df_h2$Eosinophils <= 0.75 ,]))
homo_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.75  ,]))
other_h2 <- length(rownames(df_h2[df_h2$Eosinophils ==0 &df_h2$Eosinophils <= 0.25  ,]))
doublets_h2 <- (doublets_h2*100)/all_h2
homo_h2 <- (homo_h2*100)/all_h2
other_h2 <- (other_h2*100)/all_h2
all_h2 <- (all_h2*100)/all_h2
h2_df <- data.frame(c("all","doublets","homo","other"), c(all_h2,doublets_h2,homo_h2,other_h2))
colnames(h2_df) <- c("identity","values")
h2_df$exp <- "Exp2"

all_h3 <- length(rownames(df_h3))
doublets_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.25 & df_h3$Eosinophils <= 0.75 ,]))
homo_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.75  ,]))
other_h3 <- length(rownames(df_h3[df_h3$Eosinophils ==0 &df_h3$Eosinophils <= 0.25  ,]))
doublets_h3 <- (doublets_h3*100)/all_h3
homo_h3 <- (homo_h3*100)/all_h3
other_h3 <- (other_h3*100)/all_h3
all_h3 <- (all_h3*100)/all_h3
h3_df <- data.frame(c("all","doublets","homo","other"), c(all_h3,doublets_h3,homo_h3,other_h3))
colnames(h3_df) <- c("identity","values")
h3_df$exp <- "Exp3"

all_h4 <- length(rownames(df_h4))
doublets_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.25 & df_h4$Eosinophils <= 0.75 ,]))
homo_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.75  ,]))
other_h4 <- length(rownames(df_h1[df_h4$Eosinophils ==0 &df_h4$Eosinophils <= 0.25  ,]))
doublets_h4 <- (doublets_h4*100)/all_h4
homo_h4 <- (homo_h4*100)/all_h4
other_h4 <- (other_h4*100)/all_h4
all_h4 <- (all_h4*100)/all_h4
h4_df <- data.frame(c("all","doublets","homo","other"), c(all_h4,doublets_h4,homo_h4,other_h4))
colnames(h4_df) <- c("identity","values")
h4_df$exp <- "Exp4"

all_h5 <- length(rownames(df_h5))
doublets_h5 <- length(rownames(df_h1[df_h5$Eosinophils >0.25 & df_h5$Eosinophils <= 0.75 ,]))
homo_h5 <- length(rownames(df_h5[df_h5$Eosinophils >0.75  ,]))
other_h5 <- length(rownames(df_h5[df_h5$Eosinophils ==0 &df_h5$Eosinophils <= 0.25  ,]))
doublets_h5 <- (doublets_h5*100)/all_h5
homo_h5 <- (homo_h5*100)/all_h5
other_h5 <- (other_h5*100)/all_h5
all_h5 <- (all_h5*100)/all_h5
h5_df <- data.frame(c("all","doublets","homo","other"), c(all_h5,doublets_h5,homo_h5,other_h5))
colnames(h5_df) <- c("identity","values")
h5_df$exp <- "Exp5"

all_h6 <- length(rownames(df_h6))
doublets_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.25 & df_h6$Eosinophils <= 0.75 ,]))
homo_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.75  ,]))
other_h6 <- length(rownames(df_h6[df_h6$Eosinophils ==0 &df_h6$Eosinophils <= 0.25  ,]))
doublets_h6 <- (doublets_h6*100)/all_h6
homo_h6 <- (homo_h6*100)/all_h6
other_h6 <- (other_h6*100)/all_h6
all_h6 <- (all_h6*100)/all_h6
h6_df <- data.frame(c("all","doublets","homo","other"), c(all_h6,doublets_h6,homo_h6,other_h6))
colnames(h6_df) <- c("identity","values")
h6_df$exp <- "Exp6"

all_h7 <- length(rownames(df_h7))
doublets_h7 <- length(rownames(df_h6[df_h7$Eosinophils >0.25 & df_h7$Eosinophils <= 0.75 ,]))
homo_h7 <- length(rownames(df_h7[df_h7$Eosinophils >0.75  ,]))
other_h7 <- length(rownames(df_h6[df_h7$Eosinophils ==0 &df_h7$Eosinophils <= 0.25  ,]))
doublets_h7 <- (doublets_h7*100)/all_h7
homo_h7 <- (homo_h7*100)/all_h7
other_h7 <- (other_h7*100)/all_h7
all_h7 <- (all_h7*100)/all_h7
h7_df <- data.frame(c("all","doublets","homo","other"), c(all_h7,doublets_h7,homo_h7,other_h7))
colnames(h7_df) <- c("identity","values")
h7_df$exp <- "Exp7"

df <- rbind(h1_df,h2_df)
df <- rbind(df,h3_df)
df <- rbind(df,h4_df)
df <- rbind(df,h5_df)
df <- rbind(df,h6_df)
df_scDblFinder <- rbind(df,h7_df)

### Based on scDblFinder without correction 
df <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_doublets_deconvolution_result075_no_corr.csv")

df_h1 <- df[grepl("^h1_", df$X), ]
df_h2 <- df[grepl("^h2_", df$X), ]
df_h3 <- df[grepl("^h3_", df$X), ]
df_h4 <- df[grepl("^h4_", df$X), ]
df_h5 <- df[grepl("^h5_", df$X), ]
df_h6 <- df[grepl("^h6_", df$X), ]
df_h7 <- df[grepl("^h7_", df$X), ]

all_h1 <- length(rownames(df_h1))
doublets_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.25 & df_h1$Eosinophils <= 0.75 ,]))
homo_h1 <- length(rownames(df_h1[df_h1$Eosinophils >0.75  ,]))
other_h1 <- length(rownames(df_h1[df_h1$Eosinophils ==0 &df_h1$Eosinophils <= 0.25  ,]))
doublets_h1 <- (doublets_h1*100)/all_h1
homo_h1 <- (homo_h1*100)/all_h1
other_h1 <- (other_h1*100)/all_h1
all_h1 <- (all_h1*100)/all_h1
h1_df <- data.frame(c("all","doublets","homo","other"), c(all_h1,doublets_h1,homo_h1,other_h1))
colnames(h1_df) <- c("identity","values")
h1_df$exp <- "Exp1"

all_h2 <- length(rownames(df_h2))
doublets_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.25 & df_h2$Eosinophils <= 0.75 ,]))
homo_h2 <- length(rownames(df_h2[df_h2$Eosinophils >0.75  ,]))
other_h2 <- length(rownames(df_h2[df_h2$Eosinophils ==0 &df_h2$Eosinophils <= 0.25  ,]))
doublets_h2 <- (doublets_h2*100)/all_h2
homo_h2 <- (homo_h2*100)/all_h2
other_h2 <- (other_h2*100)/all_h2
all_h2 <- (all_h2*100)/all_h2
h2_df <- data.frame(c("all","doublets","homo","other"), c(all_h2,doublets_h2,homo_h2,other_h2))
colnames(h2_df) <- c("identity","values")
h2_df$exp <- "Exp2"

all_h3 <- length(rownames(df_h3))
doublets_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.25 & df_h3$Eosinophils <= 0.75 ,]))
homo_h3 <- length(rownames(df_h3[df_h3$Eosinophils >0.75  ,]))
other_h3 <- length(rownames(df_h3[df_h3$Eosinophils ==0 &df_h3$Eosinophils <= 0.25  ,]))
doublets_h3 <- (doublets_h3*100)/all_h3
homo_h3 <- (homo_h3*100)/all_h3
other_h3 <- (other_h3*100)/all_h3
all_h3 <- (all_h3*100)/all_h3
h3_df <- data.frame(c("all","doublets","homo","other"), c(all_h3,doublets_h3,homo_h3,other_h3))
colnames(h3_df) <- c("identity","values")
h3_df$exp <- "Exp3"

all_h4 <- length(rownames(df_h4))
doublets_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.25 & df_h4$Eosinophils <= 0.75 ,]))
homo_h4 <- length(rownames(df_h4[df_h4$Eosinophils >0.75  ,]))
other_h4 <- length(rownames(df_h1[df_h4$Eosinophils ==0 &df_h4$Eosinophils <= 0.25  ,]))
doublets_h4 <- (doublets_h4*100)/all_h4
homo_h4 <- (homo_h4*100)/all_h4
other_h4 <- (other_h4*100)/all_h4
all_h4 <- (all_h4*100)/all_h4
h4_df <- data.frame(c("all","doublets","homo","other"), c(all_h4,doublets_h4,homo_h4,other_h4))
colnames(h4_df) <- c("identity","values")
h4_df$exp <- "Exp4"

all_h5 <- length(rownames(df_h5))
doublets_h5 <- length(rownames(df_h1[df_h5$Eosinophils >0.25 & df_h5$Eosinophils <= 0.75 ,]))
homo_h5 <- length(rownames(df_h5[df_h5$Eosinophils >0.75  ,]))
other_h5 <- length(rownames(df_h5[df_h5$Eosinophils ==0 &df_h5$Eosinophils <= 0.25  ,]))
doublets_h5 <- (doublets_h5*100)/all_h5
homo_h5 <- (homo_h5*100)/all_h5
other_h5 <- (other_h5*100)/all_h5
all_h5 <- (all_h5*100)/all_h5
h5_df <- data.frame(c("all","doublets","homo","other"), c(all_h5,doublets_h5,homo_h5,other_h5))
colnames(h5_df) <- c("identity","values")
h5_df$exp <- "Exp5"

all_h6 <- length(rownames(df_h6))
doublets_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.25 & df_h6$Eosinophils <= 0.75 ,]))
homo_h6 <- length(rownames(df_h6[df_h6$Eosinophils >0.75  ,]))
other_h6 <- length(rownames(df_h6[df_h6$Eosinophils ==0 &df_h6$Eosinophils <= 0.25  ,]))
doublets_h6 <- (doublets_h6*100)/all_h6
homo_h6 <- (homo_h6*100)/all_h6
other_h6 <- (other_h6*100)/all_h6
all_h6 <- (all_h6*100)/all_h6
h6_df <- data.frame(c("all","doublets","homo","other"), c(all_h6,doublets_h6,homo_h6,other_h6))
colnames(h6_df) <- c("identity","values")
h6_df$exp <- "Exp6"

all_h7 <- length(rownames(df_h7))
doublets_h7 <- length(rownames(df_h6[df_h7$Eosinophils >0.25 & df_h7$Eosinophils <= 0.75 ,]))
homo_h7 <- length(rownames(df_h7[df_h7$Eosinophils >0.75  ,]))
other_h7 <- length(rownames(df_h6[df_h7$Eosinophils ==0 &df_h7$Eosinophils <= 0.25  ,]))
doublets_h7 <- (doublets_h7*100)/all_h7
homo_h7 <- (homo_h7*100)/all_h7
other_h7 <- (other_h7*100)/all_h7
all_h7 <- (all_h7*100)/all_h7
h7_df <- data.frame(c("all","doublets","homo","other"), c(all_h7,doublets_h7,homo_h7,other_h7))
colnames(h7_df) <- c("identity","values")
h7_df$exp <- "Exp7"

df <- rbind(h1_df,h2_df)
df <- rbind(df,h3_df)
df <- rbind(df,h4_df)
df <- rbind(df,h5_df)
df <- rbind(df,h6_df)
df_scDblFinder_no_corr <- rbind(df,h7_df)

### Merge all data frames 
df_scDblFinder_no_corr$method <- "scDblFinder_no_corr"
df_scDblFinder$method <- "scDblFinder_corr"

df <- rbind(df_scDblFinder_no_corr, df_scDblFinder)
df <- df[df$identity %in% c("doublets","homo","other"),]

p <- ggplot(df, aes(x = identity, y =  values, fill = method)) + 
  geom_boxplot(outlier.shape = 16) + theme_minimal() + 
  geom_point(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = 2.5, shape = 21)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,100) +
  scale_fill_manual(values =  c( "scDblFinder_no_corr" = "#BF4BEF", "scDblFinder_corr" = "#0FD367"))
ggsave("/scratch/khandl/technical/figures/Doublet/RNA_corr_no_corr.svg", width = 10, height = 8, plot = p)
