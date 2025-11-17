######### This code generates seurat objects to assess the effect of the automatic vs. forced cell determination and mapping with and without intronic reads in BD pipeline #####   
# data from GSE282765
# Mouse AKPS CRC tumor blood and BM and healthy blood and bm 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

### automatic cell determination - intronic and exonic reads 
blood_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_blood_healthy_ST11_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

blood_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_blood_AKPS_tumor_ST12_Expression_Data.st"), 
  project = "blood_tumor", condition = "blood_tumor",3,200)

bm_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_bm_healthy_ST09_Expression_Data.st"), 
  project = "bm_wt", condition = "bm_wt",3,200)

bm_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_bm_AKPS_tumor_ST10_Expression_Data.st"), 
  project = "bm_tumor", condition = "bm_tumor",3,200)

## Merge samples
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "automatic"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Automatic_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds")

### forced cell determination - intronic and exonic reads 
blood_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_blood_healthy_ST11_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

blood_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_blood_AKPS_tumor_ST12_Expression_Data.st"), 
  project = "blood_tumor", condition = "blood_tumor",3,200)

bm_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_bm_healthy_ST09_Expression_Data.st"), 
  project = "bm_wt", condition = "bm_wt",3,200)

bm_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_bm_AKPS_tumor_ST10_Expression_Data.st"), 
  project = "bm_tumor", condition = "bm_tumor",3,200)

## Merge samples
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "forced"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_intronic_and_exonic_reads_Mm_healthy_CRC_blood_bm.rds")

### forced cell determination - exonic reads only 
blood_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_blood_healthy_ST11_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

blood_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_blood_AKPS_tumor_ST12_Expression_Data.st"), 
  project = "blood_tumor", condition = "blood_tumor",3,200)

bm_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_bm_healthy_ST09_Expression_Data.st"), 
  project = "bm_wt", condition = "bm_wt",3,200)

bm_tumor <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_blood_bm_CRC_AKPS_healthy/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_bm_AKPS_tumor_ST10_Expression_Data.st"), 
  project = "bm_tumor", condition = "bm_tumor",3,200)

## Merge samples
tumor <- merge(blood_wt, y = c(blood_tumor, bm_wt, bm_tumor),
               add.cell.ids = c( "blood_wt","blood_tumor","bm_wt","bm_tumor"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "forced"
tumor$reads <- "exonic_only"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Mm_healthy_CRC_blood_bm.rds")

