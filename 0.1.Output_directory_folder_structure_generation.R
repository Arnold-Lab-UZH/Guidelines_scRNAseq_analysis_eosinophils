########## This code creates all directories needed to save output files based on 0.config.R ##########
### Function to create the folders 
create_project_dirs <- function(base_dir = getwd()) {
  # Define the base directory 
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  
  # Top-level output roots, exactly as used in the configuration 
  data_root    <- file.path(base_dir, "Guidelines")
  results_root <- file.path(base_dir, "Guidelines", "results")
  
  # All directories defined in 0.config.R
  dirs <- c(
    file.path(data_root, "seurat_objects"),
    
    file.path(results_root, "2_data_annotation", "plots"),
    file.path(results_root, "2_data_annotation", "tables"),
    
    file.path(results_root, "3_gene_counts", "plots"),
    file.path(results_root, "3_gene_counts", "tables"),
    
    file.path(results_root, "4_dropout_and_markers", "plots"),
    file.path(results_root, "4_dropout_and_markers", "tables"),
    file.path(results_root, "4_dropout_and_markers", "plots", "Eos_marker"),
    file.path(results_root, "4_dropout_and_markers", "tables", "Eos_marker"),
    
    file.path(results_root, "5_gene_mapping", "plots"),
    file.path(results_root, "5_gene_mapping", "tables"),
    
    file.path(results_root, "6_empty_droplets", "plots"),
    file.path(results_root, "6_empty_droplets", "tables"),
    
    file.path(results_root, "7_ambient_rna", "plots"),
    file.path(results_root, "7_ambient_rna", "tables"),
    
    file.path(results_root, "8_doublet_detection", "plots"),
    file.path(results_root, "8_doublet_detection", "tables"),
    
    file.path(results_root, "9_wt_vs_il5tg", "plots"),
    file.path(results_root, "9_wt_vs_il5tg", "tables")
  )
  
  created <- vapply(dirs, function(d) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    d
  }, FUN.VALUE = character(1))
  
  message(
    length(created),
    " project output directories ensured under: ",
    base_dir
  )
  
  invisible(as.character(created))
}

### Create everything under a directory of your choice
create_project_dirs(base_dir = "/Users/handler/Downloads")

