### Function to identify doublets by the upper 99th percentile of nFeatures for each experiment separately 
annotate_doublets_based_on_upper99th_percentile <- function(
    seurat_object,
    experiment_id
){
  Idents(seurat_object) <- "experiment"
  sub <- subset(seurat_object, idents = experiment_id)
  # calculate upper 99th percentile for cutoff 
  nf <- sub$nFeature_RNA
  upper_q99 <- quantile(nf, 0.99)
  print(upper_q99)
  # extract double IDs 
  doublets <- subset(sub,subset =nFeature_RNA > upper_q99 )
  doublet_ids <- rownames(doublets@meta.data)
  
  # extract singlet IDs 
  singlets <- subset(sub,subset =nFeature_RNA <= upper_q99 )
  singlet_ids <- rownames(singlets@meta.data)
  
  sub$doublet_singlet <- NA
  sub@meta.data <- sub@meta.data %>%
    mutate(doublet_singlet = case_when(
      rownames(sub@meta.data) %in% doublet_ids ~ "doublet",
      rownames(sub@meta.data) %in% singlet_ids ~ "singlet",
      TRUE ~ NA_character_))
  return(sub)
}
