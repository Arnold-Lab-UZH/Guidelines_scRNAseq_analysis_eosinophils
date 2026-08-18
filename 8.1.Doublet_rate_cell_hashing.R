######### This code calculates the doublet rate based on sample tag cell hashing data  ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### Load Sample Tag Calls and exclude Undetermined calls 
df1 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P1_Sample_Tag_Calls.csv"),skip = 7)
df1 <- df1[!df1$Sample_Name %in% c("Undetermined"),] 
df2 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P2_Sample_Tag_Calls.csv"),skip = 7)
df2 <- df2[!df2$Sample_Name %in% c("Undetermined"),] 
df3 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P3_Sample_Tag_Calls.csv"),skip = 7)
df3 <- df3[!df3$Sample_Name %in% c("Undetermined"),] 
df4 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P4_Sample_Tag_Calls.csv"),skip = 7)
df4 <- df4[!df4$Sample_Name %in% c("Undetermined"),] 
df5 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P5_Sample_Tag_Calls.csv"),skip = 7)
df5 <- df5[!df5$Sample_Name %in% c("Undetermined"),] 
df6 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P6_Sample_Tag_Calls.csv"),skip = 7)
df6 <- df6[!df6$Sample_Name %in% c("Undetermined"),] 
df7 <- read.csv(file.path(raw_data_GSE282765_sample_tag_calls_dir,"Hs_P7_Sample_Tag_Calls.csv"),skip = 7)
df7 <- df7[!df7$Sample_Name %in% c("Undetermined"),] 

df_list <- list(df1,df2,df3,df4,df5,df6,df7)

##### Calculate multiplet rate for each dataset 
for(i in df_list) {
  all_bc <- length(rownames(i))
  multiplet_calls <- i[i$Sample_Tag %in% "Multiplet",]
  multiplet_calls <- length(rownames(multiplet_calls))
  multiplet_rate <- (multiplet_calls * 100)/all_bc
  print(multiplet_rate)
}

# experiment 7 = h6/P6, experiment 8 = H7/P7
sample <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
multiplet_rate_per_sample <- c(3.136472,6.154207,2.778757,0.5149662,0.759,4.524919,1.197831)
df <- data.frame(sample, multiplet_rate_per_sample)
df$method <- "cell_hashing"
write.csv(df,file.path(doublet_tables_dir, "cell_hashing_doublet_rate.csv"))
