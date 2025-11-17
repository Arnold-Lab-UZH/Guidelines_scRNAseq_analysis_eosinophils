######### This code generates seurat objects to assess the effect of the automatic vs. forced cell determination and mapping with and without intronic reads in BD pipeline #####   
# data from GSE182001
# Mouse il5-tg steady state bone marrow, blood, spleen, colon, small intestine, stomach 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

### automatic cell determination - intronic and exonic reads 
stomach <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_stomach_ST01_Expression_Data.st"), 
  project = "steady_state", condition = "stomach",3,200)

colon <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_colon_ST02_Expression_Data.st"), 
  project = "steady_state", condition = "colon",3,200)

small_int <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_SI_ST03_Expression_Data.st"), 
  project = "steady_state", condition = "small_int",3,200)

spleen <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_spleen_ST04_Expression_Data.st"), 
  project = "steady_state", condition = "spleen",3,200)

blood <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_blood_ST08_Expression_Data.st"), 
  project = "steady_state", condition = "blood",3,200)

bm <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_il5tg_bm_ST06_Expression_Data.st"), 
  project = "steady_state", condition = "bm",3,200)

## Merge samples
merged <- merge(stomach, y = c(colon,small_int,spleen,  blood, bm),
                add.cell.ids = c("stomach","colon","SI","spleen", "blood","bm"))
meged <- JoinLayers(merged)

## Add mitochondrial percentage per cell 
meged$percent.mt <- PercentageFeatureSet(meged, pattern = "^mt.")

meged$cell_determination <- "automatic"
meged$reads <- "intronic_and_exonic"
meged$species <- "Mm"
meged$technology <- "BD_Rhapsody"
meged$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(meged, file = "/scratch/khandl/technical/seurat_objects/Automatic_cell_determination_intronic_and_exonic_reads_Mm_il5tg_steady_state.rds")

### forced cell determination - intronic and exonic reads 
stomach <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_stomach_ST01_Expression_Data.st"), 
  project = "steady_state", condition = "stomach",3,200)

colon <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_colon_ST02_Expression_Data.st"), 
  project = "steady_state", condition = "colon",3,200)

small_int <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_SI_ST03_Expression_Data.st"), 
  project = "steady_state", condition = "small_int",3,200)

spleen <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_spleen_ST04_Expression_Data.st"), 
  project = "steady_state", condition = "spleen",3,200)

blood <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_blood_ST08_Expression_Data.st"), 
  project = "steady_state", condition = "blood",3,200)

bm <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_il5tg_bm_ST06_Expression_Data.st"), 
  project = "steady_state", condition = "bm",3,200)

## Merge samples
merged <- merge(stomach, y = c(colon,small_int,spleen,  blood, bm),
                add.cell.ids = c("stomach","colon","SI","spleen", "blood","bm"))
meged <- JoinLayers(merged)

## Add mitochondrial percentage per cell 
meged$percent.mt <- PercentageFeatureSet(meged, pattern = "^mt.")

meged$cell_determination <- "forced"
meged$reads <- "intronic_and_exonic"
meged$species <- "Mm"
meged$technology <- "BD_Rhapsody"
meged$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(meged, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_intronic_and_exonic_reads_Mm_il5tg_steady_state.rds")

### forced cell determination - exonic reads only 
stomach <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_stomach_ST01_Expression_Data.st"), 
  project = "steady_state", condition = "stomach",3,200)

colon <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_colon_ST02_Expression_Data.st"), 
  project = "steady_state", condition = "colon",3,200)

small_int <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_SI_ST03_Expression_Data.st"), 
  project = "steady_state", condition = "small_int",3,200)

spleen <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_spleen_ST04_Expression_Data.st"), 
  project = "steady_state", condition = "spleen",3,200)

blood <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_blood_ST08_Expression_Data.st"), 
  project = "steady_state", condition = "blood",3,200)

bm <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_il5tg_steady_state/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_il5tg_bm_ST06_Expression_Data.st"), 
  project = "steady_state", condition = "bm",3,200)

## Merge samples
merged <- merge(stomach, y = c(colon,small_int,spleen,  blood, bm),
                add.cell.ids = c("stomach","colon","SI","spleen", "blood","bm"))
meged <- JoinLayers(merged)

## Add mitochondrial percentage per cell 
meged$percent.mt <- PercentageFeatureSet(meged, pattern = "^mt.")

meged$cell_determination <- "forced"
meged$reads <- "exonic_only"
meged$species <- "Mm"
meged$technology <- "BD_Rhapsody"
meged$cell_enrichment  <- "Eosinophils"

## save object
saveRDS(meged, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Mm_il5tg_steady_state.rds")
