########## This code generates seurat objects to assess the effect of the automatic vs. forced cell determination and mapping with and without intronic reads in BD pipeline #####   
# data from GSE282765
# tumor and NAT from P1-P7

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

### automatic cell determination - intronic and exonic reads 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P1_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                                      "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P1_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                                        "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P2_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                      "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P2_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                                        "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P3_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                                      "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P3_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                                        "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P4_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                                      "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P4_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                                        "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P5_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                      "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P5_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                                        "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P6_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                                      "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P6_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                        "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Automatic_cell_determination_intronic_and_exonic_reads", "P7_Hs_tumor_automatic_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                      "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Automatic_cell_determination_intronic_and_exonic_reads", "P7_Hs_NAT_automatic_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                                        "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

patients$cell_determination <- "automatic"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

## save object
saveRDS(patients, file = "/scratch/khandl/technical/seurat_objects/Automatic_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds")

### Forced cell determination intronic and exonic reads 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                                      "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                                        "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                      "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                                        "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                                      "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                                        "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                                      "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                                        "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                      "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                                        "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                                      "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                                        "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                      "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                                        "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

patients$cell_determination <- "forced"
patients$reads <- "intronic_and_exonic"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

## save object
saveRDS(patients, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_and_intronic_reads_Hs_NAT_tumor.rds")

### Forced cell determination exonic reads only 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P1_Hs_tumor_forced_cell_determination_exonic_only_ST02_Expression_Data.st"), 
                                                      "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P1_Hs_NAT_forced_cell_determination_exonic_only_ST01_Expression_Data.st"), 
                                                        "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P2_Hs_tumor_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                                      "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P2_Hs_NAT_forced_cell_determination_exonic_only_ST08_Expression_Data.st"), 
                                                        "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P3_Hs_tumor_forced_cell_determination_exonic_only_ST09_Expression_Data.st"), 
                                                      "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P3_Hs_NAT_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                                                        "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P4_Hs_tumor_forced_cell_determination_exonic_only_ST04_Expression_Data.st"), 
                                                      "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P4_Hs_NAT_forced_cell_determination_exonic_only_ST03_Expression_Data.st"), 
                                                        "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P5_Hs_tumor_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                                      "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P5_Hs_NAT_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                                        "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P6_Hs_tumor_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                                      "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P6_Hs_NAT_forced_cell_determination_exonic_only_ST07_Expression_Data.st"), 
                                                        "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_tumor/Forced_cell_determination_exonic_reads_only", "P7_Hs_tumor_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                                                      "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_CRC_NAT/Forced_cell_determination_exonic_reads_only", "P7_Hs_NAT_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                                                        "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT-")

patients$cell_determination <- "forced"
patients$reads <- "exonic_only"
patients$species <- "Hs"
patients$technology <- "BD_Rhapsody"
patients$cell_enrichment  <- "CD45"

## save object
saveRDS(patients, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Hs_NAT_tumor.rds")


