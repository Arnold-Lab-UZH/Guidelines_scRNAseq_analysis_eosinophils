########## Configuration file ##########
# Source this file at the beginning of each analysis script to define
# project-relative input and output paths.

##### Define the Project base directory (directory containing this file) to load packages and function R code 
base_dir <- tryCatch({dirname(sys.frames()[[1]]$ofile)}, error = function(e) {getwd()})

##### Configure raw data directories 
scratch_dir <- Sys.getenv("scratch_dir", "/Users/handler/Downloads")
raw_data_dir      <- file.path(scratch_dir, "Guidelines/raw")

### 1.GSE282765
raw_data_GSE282765_dir      <- file.path(raw_data_dir, "1.GSE282765")

## Hs CRC tumor 
raw_data_GSE282765_Hs_CRC_tumor_dir      <- file.path(raw_data_GSE282765_dir, "Hs_CRC_tumor")
raw_data_GSE282765_Hs_CRC_tumor_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_CRC_tumor_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Hs_CRC_tumor_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Hs_CRC_tumor_dir, "Forced_cell_determination_exonic_reads_only")
raw_data_GSE282765_Hs_CRC_tumor_automatic_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_CRC_tumor_dir, "Automatic_cell_determination_intronic_and_exonic_reads")

## Hs CRC NAT 
raw_data_GSE282765_Hs_CRC_NAT_dir      <- file.path(raw_data_GSE282765_dir, "Hs_CRC_NAT")
raw_data_GSE282765_Hs_CRC_NAT_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_CRC_NAT_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Hs_CRC_NAT_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Hs_CRC_NAT_dir, "Forced_cell_determination_exonic_reads_only")
raw_data_GSE282765_Hs_CRC_NAT_automatic_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_NAT_tumor_dir, "Automatic_cell_determination_intronic_and_exonic_reads")

## Hs PB healthy 
raw_data_GSE282765_Hs_PB_healthy_dir      <- file.path(raw_data_GSE282765_dir, "Hs_PB_healthy")
raw_data_GSE282765_Hs_PB_healthy_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_PB_healthy_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Hs_PB_healthy_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Hs_PB_healthy_dir, "Forced_cell_determination_exonic_reads_only")

## Hs PB CRC
raw_data_GSE282765_Hs_PB_CRC_dir      <- file.path(raw_data_GSE282765_dir, "Hs_PB_CRC")
raw_data_GSE282765_Hs_PB_CRC_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Hs_PB_CRC_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Hs_PB_CRC_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Hs_PB_CRC_dir, "Forced_cell_determination_exonic_reads_only")

## Mm colon CRC AKPS tumor NAT disseminated
raw_data_GSE282765_Mm_colon_dir      <- file.path(raw_data_GSE282765_dir, "Mm_colon_CRC_AKPS_tumor_NAT_disseminated")
raw_data_GSE282765_Mm_colon_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Mm_colon_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Mm_colon_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Mm_colon_dir, "Forced_cell_determination_exonic_reads_only")

## Mm blood bm CRC AKPS healthy 
raw_data_GSE282765_Mm_bm_blood_dir      <- file.path(raw_data_GSE282765_dir, "Mm_blood_bm_CRC_AKPS_healthy")
raw_data_GSE282765_Mm_bm_blood_forced_intron_exon_dir      <- file.path(raw_data_GSE282765_Mm_bm_blood_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE282765_Mm_bm_blood_forced_exon_only_dir      <- file.path(raw_data_GSE282765_Mm_bm_blood_dir, "Forced_cell_determination_exonic_reads_only")
raw_data_GSE282765_Mm_bm_blood_automatic_intron_exon_dir      <- file.path(raw_data_GSE282765_Mm_bm_blood_dir, "Automatic_cell_determination_intronic_and_exonic_reads")

## Unfiltered and sample tag calls
raw_data_GSE282765_unfiltered_dir      <- file.path(raw_data_GSE282765_dir, "Unfiltered")
raw_data_GSE282765_sample_tag_calls_dir      <- file.path(raw_data_GSE282765_dir, "Sample_tag_calls")

### 2.GSE182001
raw_data_GSE182001_dir      <- file.path(raw_data_dir, "2.GSE182001")
raw_data_GSE182001_forced_intron_exon_dir      <- file.path(raw_data_GSE182001_dir, "Forced_cell_determination_intronic_and_exonic_reads")
raw_data_GSE182001_forced_exon_only_dir      <- file.path(raw_data_GSE182001_dir, "Forced_cell_determination_exonic_reads_only")

### 3.GSE256088
raw_data_GSE256088_dir      <- file.path(raw_data_dir, "3.GSE256088")

### 4.E-MTAB-14010
raw_data_E_MTAB_14010_dir      <- file.path(raw_data_dir, "4.E-MTAB-14010")

### 5.GSM7919060
raw_data_GSM7919060_dir      <- file.path(raw_data_dir, "5.GSM7919060")

### 6.GSE276583
raw_data_GSE276583_dir      <- file.path(raw_data_dir, "6.GSE276583")

### 7.GSE216189
raw_data_GSE216189_dir      <- file.path(raw_data_dir, "7.GSE216189")

### 8.GSE175930
raw_data_GSE175930_dir      <- file.path(raw_data_dir, "8.GSE175930")

##### Configure data output directories (relative to repository)
data_dir          <- Sys.getenv("data_dir", "/data/khandl")
seurat_objects_dir <- file.path(data_dir, "Guidelines/seurat_objects")

##### Configure results output directories (plots and tables)
results_dir <- file.path(data_dir, "Guidelines/results")

### 2 Data annotation  
annotation_plots_dir  <- file.path(results_dir, "2_data_annotation", "plots")
annotation_tables_dir <- file.path(results_dir, "2_data_annotation", "tables")

### 3 Gene counts per cell type
gene_counts_plots_dir  <- file.path(results_dir, "3_gene_counts", "plots")
gene_counts_tables_dir <- file.path(results_dir, "3_gene_counts", "tables")

# 4. Gene dropout rates and marker gene expression in eosinophils
dropout_plots_dir  <- file.path(results_dir, "4_dropout_and_markers", "plots")
dropout_tables_dir <- file.path(results_dir, "4_dropout_and_markers", "tables")
eos_marker_plots_dir <- file.path(dropout_plots_dir, "Eos_marker")
eos_marker_tables_dir <- file.path(dropout_tables_dir, "Eos_marker")

# 5. Gene mapping to intronic and exonic reads or exonic reads only
gene_mapping_plots_dir  <- file.path(results_dir, "5_gene_mapping", "plots")
gene_mapping_tables_dir <- file.path(results_dir, "5_gene_mapping", "tables")

# 6. Empty wells detection
empty_droplets_plots_dir  <- file.path(results_dir, "6_empty_droplets", "plots")
empty_droplets_tables_dir <- file.path(results_dir, "6_empty_droplets", "tables")

# 7. Ambient RNA detection and correction
ambient_rna_plots_dir  <- file.path(results_dir, "7_ambient_rna", "plots")
ambient_rna_tables_dir <- file.path(results_dir, "7_ambient_rna", "tables")

# 8. Doublet detection
doublet_plots_dir  <- file.path(results_dir, "8_doublet_detection", "plots")
doublet_tables_dir <- file.path(results_dir, "8_doublet_detection", "tables")

# 9. Comparison of eosinophils from Il5-tg and wild type BM and colon
wt_il5tg_plots_dir  <- file.path(results_dir, "9_wt_vs_il5tg", "plots")
wt_il5tg_tables_dir <- file.path(results_dir, "9_wt_vs_il5tg", "tables")
