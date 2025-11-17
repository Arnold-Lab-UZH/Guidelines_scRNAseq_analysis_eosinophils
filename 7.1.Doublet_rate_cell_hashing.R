######### This code compares the mesians of number of features genes of annotated cell types across datasets  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### read R objects 
df1 <- read.csv("/scratch/khandl/Sample_tag_calls/Mm_colon_healthy_wt_phil_Sample_Tag_Calls.csv",skip = 7)
df1 <- df1[!df1$Sample_Name %in% c("Undetermined"),] 
df2 <- read.csv("/scratch/khandl/Sample_tag_calls/Mm_tumor_phil_dissemianted_phil_wt_NAT_wt_Sample_Tag_Calls.csv",skip = 7)
df2 <- df2[!df2$Sample_Name %in% c("Undetermined"),] 
df3 <- read.csv("/scratch/khandl/Technical3/Mm_BM_blood_tumor_healthy_Sample_Tag_Calls.csv",skip = 7)
df3 <- df3[!df3$Sample_Name %in% c("Undetermined"),] 
df4 <- read.csv("/scratch/khandl/Technical3/Mm_il5tg_bm_blood_Sample_Tag_Calls.csv",skip = 7)
df4 <- df4[!df4$Sample_Name %in% c("Undetermined"),] 
df5 <- read.csv("/scratch/khandl/Technical3/Mm_il5tg_stomach_colon_small_int_spleen_Sample_Tag_Calls.csv",skip = 7)
df5 <- df5[!df5$Sample_Name %in% c("Undetermined"),] 
df6 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp1_Sample_Tag_Calls.csv",skip = 7)
df6 <- df6[!df6$Sample_Name %in% c("Undetermined"),] 
df7 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp2_Sample_Tag_Calls.csv",skip = 7)
df7 <- df7[!df7$Sample_Name %in% c("Undetermined"),] 
df8 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp3_Sample_Tag_Calls.csv",skip = 7)
df8 <- df8[!df8$Sample_Name %in% c("Undetermined"),] 
df9 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp4_Sample_Tag_Calls.csv",skip = 7)
df9 <- df9[!df9$Sample_Name %in% c("Undetermined"),] 
df10 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp5_Sample_Tag_Calls.csv",skip = 7)
df110 <- df10[!df10$Sample_Name %in% c("Undetermined"),] 
df11 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp6_Sample_Tag_Calls.csv",skip = 7)
df11 <- df11[!df11$Sample_Name %in% c("Undetermined"),] 
df12 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp7_Sample_Tag_Calls.csv",skip = 7)
df12 <- df12[!df12$Sample_Name %in% c("Undetermined"),] 
df13 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp_PB_healthy_Sample_Tag_Calls.csv",skip = 7)
df13 <- df13[!df13$Sample_Name %in% c("Undetermined"),] 

df_list <- list(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13)

##### calculate multiplet reate for each dataset 
for(i in df_list) {
  all_bc <- length(rownames(i))
  multiplet_calls <- i[i$Sample_Tag %in% "Multiplet",]
  multiplet_calls <- length(rownames(multiplet_calls))
  multiplet_rate <- (multiplet_calls * 100)/all_bc
  print(multiplet_rate)
}

##### Doublet rate for Human Exp 1-7
df6 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp1_Sample_Tag_Calls.csv",skip = 7)
df6 <- df6[!df6$Sample_Name %in% c("Undetermined"),] 
df7 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp2_Sample_Tag_Calls.csv",skip = 7)
df7 <- df7[!df7$Sample_Name %in% c("Undetermined"),] 
df8 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp3_Sample_Tag_Calls.csv",skip = 7)
df8 <- df8[!df8$Sample_Name %in% c("Undetermined"),] 
df9 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp4_Sample_Tag_Calls.csv",skip = 7)
df9 <- df9[!df9$Sample_Name %in% c("Undetermined"),] 
df10 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp5_Sample_Tag_Calls.csv",skip = 7)
df110 <- df10[!df10$Sample_Name %in% c("Undetermined"),] 
df11 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp6_Sample_Tag_Calls.csv",skip = 7)
df11 <- df11[!df11$Sample_Name %in% c("Undetermined"),] 
df12 <- read.csv("/scratch/khandl/Sample_tag_calls/Hs_exp7_Sample_Tag_Calls.csv",skip = 7)
df12 <- df12[!df12$Sample_Name %in% c("Undetermined"),] 

df_list <- list(df6,df7,df8,df9,df10,df11,df12)

##### calculate multiplet reate for each dataset 
for(i in df_list) {
  all_bc <- length(rownames(i))
  multiplet_calls <- i[i$Sample_Tag %in% "Multiplet",]
  multiplet_calls <- length(rownames(multiplet_calls))
  multiplet_rate <- (multiplet_calls * 100)/all_bc
  print(multiplet_rate)
}

sample <- c("Exp1","Exp2","Exp3","Exp4","Exp5","Exp7","Exp8")
multiplet_rate_per_sample <- c(3.136472,6.154207,2.778757,0.5149662,0.759,4.524919,1.197831)
df <- data.frame(sample, multiplet_rate_per_sample)
df$method <- "cell_hashing"
write.csv(df,"/scratch/khandl/technical/figures/Doublet/cell_hashing_doublet_rate.csv")
