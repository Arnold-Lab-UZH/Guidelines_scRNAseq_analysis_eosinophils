### Function to read in data generated from BD Seven Bridges platform 
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Bioproduct, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Bioproduct")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}

data_to_sparse_matrix_unfiltered <- function(data.st_file_path) {
  # 1. Read the BD Rhapsody .st file
  # Expect columns: Cell_Index, Bioproduct, RSEC_Adjusted_Molecules
  input <- read.table(data.st_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 2. Clean and prepare
  input <- input %>%
    filter(!is.na(Cell_Index), !is.na(Bioproduct)) %>%
    mutate(
      Cell_Index = as.character(Cell_Index),
      Bioproduct = as.character(Bioproduct)
    )
  
  # 3. Convert to sparse matrix directly
  # factor() assigns integer indices for genes (i) and cells (j)
  gene_levels <- unique(input$Bioproduct)
  cell_levels <- unique(input$Cell_Index)
  
  sparse_mat <- sparseMatrix(
    i = as.integer(factor(input$Bioproduct, levels = gene_levels)),
    j = as.integer(factor(input$Cell_Index, levels = cell_levels)),
    x = input$RSEC_Adjusted_Molecules,
    dimnames = list(gene_levels, cell_levels)
  )
  return(sparse_mat)
}

### Seurat object generation for Homo sapiens data
create_seurat_Hs_data <- function(
    path_to_st_file,
    project,
    min.cells,
    min.features,
    condition_oi,
    tissue_oi,
    experiment_oi,
    phenotype_oi
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  seurat_object <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  seurat_object$condition <- condition_oi
  seurat_object$tissue <- tissue_oi
  seurat_object$experiment <- experiment_oi
  seurat_object$phenotype <- phenotype_oi
  return(seurat_object)
}

create_seurat_Hs_data_from_sparse_matrix <- function(
    input_matrix,
    project,
    min.cells,
    min.features,
    condition_oi,
    tissue_oi,
    experiment_oi,
    phenotype_oi
) {
  seurat_object <-CreateSeuratObject(input_matrix, 
                                     project = project,
                                     min.cells = min.cells,
                                     min.features = min.features)
  seurat_object$condition <- condition_oi
  seurat_object$tissue <- tissue_oi
  seurat_object$experiment <- experiment_oi
  seurat_object$phenotype <- phenotype_oi
  return(seurat_object)
}

### Seurat object generation for Mus musculus data 
create_seurat_Mm_data <- function(
    path_to_st_file,
    project,
    condition, 
    min.cells,
    min.features
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  condition_sample$condition <- condition
  return(condition_sample)
}

create_seurat_Mm_data_from_sparse_matrix <- function(
    input_matrix,
    project,
    condition, 
    min.cells,
    min.features
) {
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  condition_sample$condition <- condition
  return(condition_sample)
}

### Seurat object generation for 10X data 
create_seurat_10X_structured_data <- function(
    path_to_data_file,
    project,
    min.cells,
    min.features,
    condition_oi,
    tissue_oi,
    experiment_oi,
    phenotype_oi
) {
  dat <- Read10X(data.dir = path_to_data_file)
  dat <- dat[!is.na(rownames(dat)), ]
  seurat_object <-CreateSeuratObject(dat,  project = project,min.cells = min.cells, min.features = min.features)
  seurat_object$condition <- condition_oi
  seurat_object$tissue <- tissue_oi
  seurat_object$experiment <- experiment_oi
  seurat_object$phenotype <- phenotype_oi
  return(seurat_object)
}

### Seurat object generation from zUMI outputs 
create_seurat_zUMI_outputs <- function(
    count_matrix,
    project,
    min.cells,
    min.features,
    condition, 
    batch,
    tissue_oi,
    phenotype_oi
) {
  zUMI_output <- readRDS(count_matrix)
  #use all exon UMI counts from zUMI output 
  dataframe<-as.data.frame(as.matrix(zUMI_output$umicount$exon$all))%>%
    as.matrix(.)
  dataframe <- dataframe[2:(length(rownames(dataframe))),]
  #add column name "ensemble_gene_id_version" to the column of Ensemble IDs in zUMI output count matrix 
  dataframe<-mutate(as.data.frame(dataframe),ensembl_gene_id_version=rownames(dataframe))
  #compare ensemble IDs of zUMI count matrix with gene ID data frame to match ensemble IDs and to rename them to external gene names 
  gene_ids <- gene_ids_mouse 
  gene_ids$X <- NULL
  
  join<-dataframe%>%
    left_join(dplyr::select(gene_ids,1:2))
  #remove duplicates 
  join<-join[!duplicated(join$external_gene_name),]
  join[is.na(join)]<-0 #make all empty value to zero
  #add external gene names as rownames (instead of ensemble IDs before)
  rownames(join)<-join$external_gene_name
  #remove unnessesary columns 
  join<-dplyr::select(join,-ensembl_gene_id_version,-external_gene_name)
  seuratObject <-CreateSeuratObject(join, 
                                    project = project,
                                    min.cells = min.cells,
                                    min.features = min.features)
  seuratObject$condition <- condition
  seuratObject$batch <- batch 
  seuratObject$tissue <- tissue_oi
  seuratObject$phenotype <- phenotype_oi
  return(seuratObject)
}


