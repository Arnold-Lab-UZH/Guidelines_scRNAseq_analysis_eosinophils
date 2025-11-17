######### This code generates seurat objects to assess the effect of the automatic vs. forced cell determination and mapping with and without intronic reads in BD pipeline #####   
# data from GSE282765
# Mouse AKPS CRC tumor, NAT, disseminated and colon healthy 

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_preprocessing.R")

### automatic cell determination - intronic and exonic reads 
tumor_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_AKPS_tumor_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_AKPS_disseminated_ST07_Expression_Data.st"), 
  project = "disseminated_wt", condition = "disseminated_wt",3,200)

adult_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_colon_ST04_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

adjacent_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Automatic_cell_determination_intronic_and_exonic_reads", "Automatic_cell_determination_intronic_and_exonic_Mm_AKPS_NAT_ST03_Expression_Data.st"), 
  project = "adjacent_colon_wt", condition = "adjacent_colon_wt",3,200)

## Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt,adjacent_colon_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt","adjacent_colon_wt"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "automatic"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Automatic_cell_determination_intronic_and_exonic_reads_Mm_CRC_AKPS_tumor_NAT_diss.rds")

### forced cell determination - intronic and exonic reads 
tumor_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_tumor_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_disseminated_ST07_Expression_Data.st"), 
  project = "disseminated_wt", condition = "disseminated_wt",3,200)

adult_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_colon_ST04_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

adjacent_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_intronic_and_exonic_reads", "Forced_cell_determination_intronic_and_exonic_Mm_AKPS_NAT_ST03_Expression_Data.st"), 
  project = "adjacent_colon_wt", condition = "adjacent_colon_wt",3,200)

## Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt,adjacent_colon_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt","adjacent_colon_wt"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "forced"
tumor$reads <- "intronic_and_exonic"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_intronic_and_exonic_reads_Mm_CRC_AKPS_tumor_NAT_diss.rds")

### forced cell determination - exonic reads only 
tumor_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical", "Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_tumor_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_disseminated_ST07_Expression_Data.st"), 
  project = "disseminated_wt", condition = "disseminated_wt",3,200)

adult_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_colon_ST04_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

adjacent_colon_wt <- create_seurat_Mm_data(
  path_to_st_file = file.path("/scratch/khandl/technical","Mm_colon_CRC_AKPS_tumor_NAT_disseminated/Forced_cell_determination_exonic_reads_only", "Forced_cell_determination_exonic_only_Mm_AKPS_NAT_ST03_Expression_Data.st"), 
  project = "adjacent_colon_wt", condition = "adjacent_colon_wt",3,200)

## Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt,adjacent_colon_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt","adjacent_colon_wt"))
tumor <- JoinLayers(tumor)

## Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

tumor$cell_determination <- "forced"
tumor$reads <- "exonic_only"
tumor$species <- "Mm"
tumor$technology <- "BD_Rhapsody"
tumor$cell_enrichment  <- "CD45"

## save object
saveRDS(tumor, file = "/scratch/khandl/technical/seurat_objects/Forced_cell_determination_exonic_reads_only_Mm_CRC_AKPS_tumor_NAT_diss.rds")



