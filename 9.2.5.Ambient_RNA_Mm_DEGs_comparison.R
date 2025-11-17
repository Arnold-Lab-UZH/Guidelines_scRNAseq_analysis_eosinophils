########## This code quantifies ambient RNA content  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.3.Functions_gene_expression.R")

# GSE282765 and GSE182001
# take forced cell determination and intronic + exonic reads

##### read R objects 
obj_RNA_scCDC <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_scCDC.rds")
obj_decontX <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_decontX.rds")
obj_SoupX <- readRDS( "/scratch/khandl/technical/seurat_objects/Mm_tumor_colon_NAT_diss_forced_cell_determination_with_intronic_reads_annotated_SoupX.rds")

Idents(obj_RNA_scCDC) <- "annotation"
obj_RNA_scCDC <- subset(obj_RNA_scCDC, idents = c("Eosinophils"))
Idents(obj_decontX) <- "annotation"
obj_decontX <- subset(obj_decontX, idents = c("Eosinophils"))
Idents(obj_SoupX) <- "annotation"
obj_SoupX <- subset(obj_SoupX, idents = c("Eosinophils"))

obj_RNA_scCDC <- JoinLayers(obj_RNA_scCDC, assay = "RNA")
obj_RNA_scCDC <- JoinLayers(obj_RNA_scCDC, assay = "Corrected")
obj_decontX <- JoinLayers(obj_decontX)
obj_SoupX <- JoinLayers(obj_SoupX)

##### DEGs between eosinopohils from tumor and healhty colon across different decontamination assays 
Idents(obj_RNA_scCDC) <- "condition"
DEG_to_csv_two_cond(obj_RNA_scCDC,"RNA", "tumor_wt","adult_colon_wt",FALSE,0.25,"/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_RNA.csv")
DEG_to_csv_two_cond(obj_RNA_scCDC,"Corrected", "tumor_wt","adult_colon_wt",FALSE,0.25,"/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_scCDC.csv")

Idents(obj_SoupX) <- "condition"
DEG_to_csv_two_cond(obj_SoupX,"RNA", "tumor_wt","adult_colon_wt",FALSE,0.25,"/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_SoupX.csv")

Idents(obj_decontX) <- "condition"
DEG_to_csv_two_cond(obj_decontX,"RNA", "tumor_wt","adult_colon_wt",FALSE,0.25,"/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_decontX.csv")

##### extract signfiicant genes high in tumor or colon and then intersect them 
df_RNA <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_RNA.csv")
df_RNA <- df_RNA[df_RNA$p_val_adj <= 0.05,]

df_decontX <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_decontX.csv")
df_decontX <- df_decontX[df_decontX$p_val_adj <= 0.05,]

df_SoupX <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_SoupX.csv")
df_SoupX <- df_SoupX[df_SoupX$p_val_adj <= 0.05,]

df_scCDC <- read.csv("/scratch/khandl/technical/figures/Ambient_RNA/DEGs_eos_tumor_vs_colon_healthy_assay_scCDC.csv")
df_scCDC <- df_scCDC[df_scCDC$p_val_adj <= 0.05,]

x <- list( "RNA" = df_RNA$X, "decontX" =df_decontX$X,
           "SoupX" = df_SoupX$X,"scCDC" = df_scCDC$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#EABBC4")

## comparison RNA vs. decontX 
x <- list( "RNA" = df_RNA$X, "decontX" =df_decontX$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#159C48")

shared <- intersect(df_RNA$X, df_decontX$X)
shared
RNA_unique <- df_RNA$X[!df_RNA$X %in% shared]
RNA_unique
decontX_unique <- df_decontX$X[!df_decontX$X %in% shared] 
decontX_unique

## comparison RNA vs. SoupX 
x <- list( "RNA" = df_RNA$X, "SoupX" =df_SoupX$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#8D3B96")

shared <- intersect(df_RNA$X, df_SoupX$X)
shared
RNA_unique <- df_RNA$X[!df_RNA$X %in% shared]
RNA_unique
SoupX_unique <- df_SoupX$X[!df_SoupX$X %in% shared] 
SoupX_unique

## comparison RNA vs. scCDC 
x <- list( "RNA" = df_RNA$X, "scCDC" =df_scCDC$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#21A1DB")

shared <- intersect(df_RNA$X, df_scCDC$X)
shared
RNA_unique <- df_RNA$X[!df_RNA$X %in% shared]
RNA_unique
scCDC_unique <- df_scCDC$X[!df_scCDC$X %in% shared] 
scCDC_unique

