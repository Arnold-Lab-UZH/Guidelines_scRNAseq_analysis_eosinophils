########## This ccode extracts barcodes with lower than 200 genes and larger 25% mito and adds it to the annotated forced with introns object ##########
## on the example of Mm tumor, diss, colon, NAT  GSE282765

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

##### load annoated object 

###lower 200 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P1",3,0,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,0, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,0,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,0,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,0,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,0,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,0,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,0,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,0,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,0,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,0,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,0,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,0,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,0,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
lower200 <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
lower200 <- JoinLayers(lower200)
## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
lower200$percent.mt <- PercentageFeatureSet(lower200, pattern = "^MT-")
lower200 <- subset(lower200, subset = nFeature_RNA < 200)
lower200$annotation <- "features_lower200"

### larger 25 mito 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
larger25 <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
larger25$percent.mt <- PercentageFeatureSet(larger25, pattern = "^MT-")
larger25 <- subset(larger25, subset = percent.mt >= 25)
larger25$annotation <- "mito_larger25"

##### merge larger25 and lower200
tumor <- merge(lower200, y = c(larger25),
               add.cell.ids = c("1","2"))
tumor <- JoinLayers(tumor)
tumor <- NormalizeData(tumor)

### save object 
saveRDS( tumor, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_lower200features_larger25mito.rds")

