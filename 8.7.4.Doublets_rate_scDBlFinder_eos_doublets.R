######### This code analyzes eosinophil contributing doublets based on position in UMAP and marker gene expression from scDBlFinder output ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load seurat objects 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Singlets_and_Doublets_Hs_NAT_tumor_PB_P1_to_P7_scCDC_annotated.rds")

##### Add doublet/singlet condition
scDblFinder_df <- read.csv("/scratch/khandl/technical/figures/Doublet/scDblFinder_wo_doublets_deconvolution_result075_RNA_corr.csv")
doublet_cell_ids <- scDblFinder_df$X
all_bcs <- rownames(obj@meta.data)
singlet_cell_ids <- all_bcs[!all_bcs %in% doublet_cell_ids]
eos_doublets_hetero <- (scDblFinder_df[scDblFinder_df$Eosinophils >0.25 & scDblFinder_df$Eosinophils <= 0.75,])$X
eos_doublets_homo <- (scDblFinder_df[scDblFinder_df$Eosinophils >0.75 ,])$X
doublet_other <- doublet_cell_ids[!doublet_cell_ids %in% c(eos_doublets_hetero,eos_doublets_homo)]

obj$doublet_cond <- NA
obj@meta.data <- obj@meta.data %>%
  mutate(doublet_cond = case_when(
    rownames(obj@meta.data) %in% singlet_cell_ids ~ "Singlet",
    rownames(obj@meta.data)  %in% doublet_other ~ "Doublet",
    rownames(obj@meta.data)  %in% eos_doublets_hetero ~ "Eos_hetero",
    rownames(obj@meta.data)  %in% eos_doublets_homo ~ "Eos_homo",
    TRUE ~ NA_character_))
table(obj$doublet_cond)
DimPlot(obj, group.by = "doublet_cond", label = TRUE)

### Highlight hetero and homotypic doublets 
Idents(obj) <- "doublet_cond"
eos_homo <- WhichCells(obj, idents = c("Eos_homo"))
eos_hetero <- WhichCells(obj, idents = c("Eos_hetero"))
p <- DimPlot(obj, label=T, group.by="doublet_cond", cells.highlight= list(eos_homo, eos_hetero), 
        cols.highlight = c( "#270CEF","#EFA40F"), cols= "#A39F9F",raster = FALSE)
ggsave("/scratch/khandl/technical/figures/Doublet/UMAP_hetero_homo_doublets_eos.svg", width = 8, height = 8, plot = p)

##### Marker gene expression in eos doublets 
Idents(obj) <- "doublet_cond"
sub <- subset(obj, idents = c("Eos_hetero","Eos_homo","Doublet"))

p <- DotPlot(sub, features = c("CLC","CCR3","ALOX15","ADGRE1","DACH1", "CD3E","ICOS","CD4", "C1QC","C1QB","VCAN","FN1","FCGR3B"),dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/technical/figures/Doublet/Eos_scDblFinder_homo_hetero.svg", width = 10, height = 6, plot = p)
