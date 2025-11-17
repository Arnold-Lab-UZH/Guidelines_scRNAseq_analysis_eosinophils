######### This code generates seurat objects to assess the effect of the automatic vs. forced cell determination and mapping with and without intronic reads in BD pipeline #####   
# data from GSE282765
# peripheral blood healthy and CRC patient 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

#scp -r /Volumes/LaCie_Drive/technical khandl@cluster.s3it.uzh.ch:scratch/ 

### automatic cell determination - intronic and exonic reads 
## healthy individuals 
H1 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H1_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                                "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H2_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                                "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H3_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                                "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H4_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                                "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H5_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                                "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H6_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")
H7 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H7_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                                "H7",3,200,  "H7_blood","blood_healthy","Exp1","healthy")
H8 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H8_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                                "H8",3,200,  "H8_blood","blood_healthy","Exp2","healthy")
H9 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "H9_Hs_PB_healthy_automatic_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                "H9",3,200,  "H9_blood","blood_healthy","Exp3","healthy")

## patient blood 
P1_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P1_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST3_Expression_Data.st"), 
                                                      "P1",3,200,  "P1_blood","blood_patient","Exp1","patient")
P2_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P2_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                                      "P2",3,200,  "P2_blood","blood_patient","Exp2","patient")
P3_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P3_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                                      "P3",3,200,  "P3_blood","blood_patient","Exp3","patient")
P4_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P4_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                                      "P4",3,200,  "P4_blood","blood_patient","Exp4","patient")
P5_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P5_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                                      "P5",3,200,  "P5_blood","blood_patient","Exp5","patient")
P6_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P6_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                                      "P6",3,200,  "P6_blood","blood_patient","Exp7","patient")
P7_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Automatic_cell_determination_intronic_and_exonic_reads", "P7_Hs_PB_automatic_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                                      "P7",3,200,  "P7_blood","blood_patient","Exp8","patient")

## Merge samples
blood <- merge(H1, y = c(H2,H3,H4,H5,H6,H7,H8,H9,P1_blood,P2_blood,P3_blood,P4_blood,P5_blood,P6_blood,P7_blood),
                       add.cell.ids = c("b1", "b2","b3","b4", "b5","b6","b7","b8","b9","pb1", "pb2","pb3","pb4", "pb5","pb6","pb7"))
blood <- JoinLayers(blood)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

blood$cell_determination <- "automatic"
blood$reads <- "intronic_and_exonic"
blood$species <- "Hs"
blood$technology <- "BD_Rhapsody"
blood$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(blood, file = "/scratch/khandl/technical/seurat_objects/Automatic_cell_determination_exonic_and_intronic_reads_Hs_PB_healthy_and_CRC.rds")

### Forced cell determination intronic and exonic reads 
## healthy individuals 
H1 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H1_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                                "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H2_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                                "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H3_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                                "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H4_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                                "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H5_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                                "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H6_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")
H7 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H7_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                                "H7",3,200,  "H7_blood","blood_healthy","Exp1","healthy")
H8 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H8_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                                "H8",3,200,  "H8_blood","blood_healthy","Exp2","healthy")
H9 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_intronic_and_exonic_reads", "H9_Hs_PB_healthy_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                                "H9",3,200,  "H9_blood","blood_healthy","Exp3","healthy")

## patient blood 
P1_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P1_Hs_PB_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                  "P1",3,200,  "P1_blood","blood_patient","Exp1","patient")
P2_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P2_Hs_PB_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P2",3,200,  "P2_blood","blood_patient","Exp2","patient")
P3_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P3_Hs_PB_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                  "P3",3,200,  "P3_blood","blood_patient","Exp3","patient")
P4_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P4_Hs_PB_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P4",3,200,  "P4_blood","blood_patient","Exp4","patient")
P5_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P5_Hs_PB_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                  "P5",3,200,  "P5_blood","blood_patient","Exp5","patient")
P6_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P6_Hs_PB_forced_cell_determination_intronic_and_exonic_ST05_Expression_Data.st"), 
                                  "P6",3,200,  "P6_blood","blood_patient","Exp7","patient")
P7_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_intronic_and_exonic_reads", "P7_Hs_PB_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                  "P7",3,200,  "P7_blood","blood_patient","Exp8","patient")

## Merge samples
blood <- merge(H1, y = c(H2,H3,H4,H5,H6,H7,H8,H9,P1_blood,P2_blood,P3_blood,P4_blood,P5_blood,P6_blood,P7_blood),
               add.cell.ids = c("b1", "b2","b3","b4", "b5","b6","b7","b8","b9","pb1", "pb2","pb3","pb4", "pb5","pb6","pb7"))
blood <- JoinLayers(blood)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

blood$cell_determination <- "forced"
blood$reads <- "intronic_and_exonic"
blood$species <- "Hs"
blood$technology <- "BD_Rhapsody"
blood$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(blood, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_and_intronic_reads_Hs_PB_healthy_and_CRC.rds")

### Forced cell determination exonic reads only 
## healthy individuals 
H1 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H1_Hs_PB_healthy_forced_cell_determination_exonic_only_ST01_Expression_Data.st"), 
                                                "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H2_Hs_PB_healthy_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                                                "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H3_Hs_PB_healthy_forced_cell_determination_exonic_only_ST09_Expression_Data.st"), 
                                                "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H4_Hs_PB_healthy_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                                                "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H5_Hs_PB_healthy_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                                                "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H6_Hs_PB_healthy_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                                                "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")
H7 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H7_Hs_PB_healthy_forced_cell_determination_exonic_only_ST04_Expression_Data.st"), 
                                                "H7",3,200,  "H7_blood","blood_healthy","Exp1","healthy")
H8 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H8_Hs_PB_healthy_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                                                "H8",3,200,  "H8_blood","blood_healthy","Exp2","healthy")
H9 <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_healthy/Forced_cell_determination_exonic_reads_only", "H9_Hs_PB_healthy_forced_cell_determination_exonic_only_ST12_Expression_Data.st"), 
                                                "H9",3,200,  "H9_blood","blood_healthy","Exp3","healthy")

## patient blood 
P1_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P1_Hs_PB_forced_cell_determination_exonic_only_ST03_Expression_Data.st"), 
                                  "P1",3,200,  "P1_blood","blood_patient","Exp1","patient")
P2_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P2_Hs_PB_forced_cell_determination_exonic_only_ST06_Expression_Data.st"), 
                                  "P2",3,200,  "P2_blood","blood_patient","Exp2","patient")
P3_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P3_Hs_PB_forced_cell_determination_exonic_only_ST11_Expression_Data.st"), 
                                  "P3",3,200,  "P3_blood","blood_patient","Exp3","patient")
P4_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P4_Hs_PB_forced_cell_determination_exonic_only_ST02_Expression_Data.st"), 
                                  "P4",3,200,  "P4_blood","blood_patient","Exp4","patient")
P5_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P5_Hs_PB_forced_cell_determination_exonic_only_ST08_Expression_Data.st"), 
                                  "P5",3,200,  "P5_blood","blood_patient","Exp5","patient")
P6_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P6_Hs_PB_forced_cell_determination_exonic_only_ST05_Expression_Data.st"), 
                                  "P6",3,200,  "P6_blood","blood_patient","Exp7","patient")
P7_blood <- create_seurat_Hs_data(file.path("/scratch/khandl/technical", "Hs_PB_CRC/Forced_cell_determination_exonic_reads_only", "P7_Hs_PB_forced_cell_determination_exonic_only_ST10_Expression_Data.st"), 
                                  "P7",3,200,  "P7_blood","blood_patient","Exp8","patient")

## Merge samples
blood <- merge(H1, y = c(H2,H3,H4,H5,H6,H7,H8,H9,P1_blood,P2_blood,P3_blood,P4_blood,P5_blood,P6_blood,P7_blood),
               add.cell.ids = c("b1", "b2","b3","b4", "b5","b6","b7","b8","b9","pb1", "pb2","pb3","pb4", "pb5","pb6","pb7"))
blood <- JoinLayers(blood)

## Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
blood$percent.mt <- PercentageFeatureSet(blood, pattern = "^MT-")

blood$cell_determination <- "forced"
blood$reads <- "exonic_only"
blood$species <- "Hs"
blood$technology <- "BD_Rhapsody"
blood$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(blood, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Hs_PB_healthy_and_CRC.rds")



