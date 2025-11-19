##### Functions to compare cell type proportions 
### create table 
create_table_cell_type_prop <- function(
    seurat_object, 
    ident_of_interest1,
    ident_of_interest2,
    output_file_path,
    sample_name
) {
  ann_tab<-  table(seurat_object@meta.data[[ident_of_interest1]], seurat_object@meta.data[[ident_of_interest2]])
  ann_tab <- cbind(ann_tab, Total = rowSums(ann_tab))
  ann_tab <- as.data.frame(ann_tab)
  ann_tab_pct = lapply(ann_tab[,], function(x) {
    x/ann_tab$Total})
  ann_tab_pct <- as.data.frame(ann_tab_pct)
  rownames(ann_tab_pct) <- rownames(ann_tab)
  write.csv(ann_tab_pct, file = paste0(output_file_path,sample_name,"_proportions_",ident_of_interest1,"_",ident_of_interest2,".csv"))
}

### change the table to plot in barplot 
create_table_cell_type_prop_for_plot_one_sample <- function(
    df_cell_type_prop_per_sample, 
    sample_name,
    n_cell_types_dim
){
  df <- df_cell_type_prop_per_sample[df_cell_type_prop_per_sample$X %in% sample_name,]
  df <- t(df)
  df <- as.data.frame(df)
  colnames(df) <- NULL
  df <- df[n_cell_types_dim,]
  df$cell_types <- rownames(df)
  colnames(df) <- c("proportion","cell_types")
  rownames(df) <- NULL
  df$sample <- sample_name
  return(df)
}

create_table_cell_type_prop_table_for_plot <- function(
    df_cell_type_prop_per_sample, 
    n_cell_types_dim
) {
  df_cell_type_prop_per_sample$Total <- NULL
  
  sample_names <- df_cell_type_prop_per_sample$X
  df <- NULL
  
  for (i in sample_names){
    df1 <- create_table_cell_type_prop_for_plot_one_sample(df_cell_type_prop_per_sample,i,n_cell_types_dim) 
    df <- rbind(df, df1)
  }
  
  df$proportion <- as.numeric(df$proportion)
  return(df)
}
