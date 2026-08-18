########### This code analyses DEGs between il5-tg and wt colon and BM eosinophils ##########
### Datasets used: GSE182001, GSE282765

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.3.Functions_gene_expression.R"))

##### Load Seurat object
obj <- readRDS( file.path(seurat_objects_dir,"wt_il5tg_bm_colon_annotated.rds"))

obj <- JoinLayers(obj)

##### Per tissue 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
Idents(sub) <- "genotype"
DEG_to_csv_two_cond(sub,"RNA", "il5tg","wt",FALSE,0.25,file.path(wt_il5tg_tables_dir, "DEGs_colon_il5tg_vs_wt.csv"))


sub <- subset(obj, idents = "bm")
Idents(sub) <- "genotype"
DEG_to_csv_two_cond(sub,"RNA", "il5tg","wt",FALSE,0.25,file.path(wt_il5tg_tables_dir, "DEGs_bm_il5tg_vs_wt.csv"))

obj$tissue_cond <- paste0(obj$tissue, "_",obj$condition)

##### Plot average expression of GOI between tissue and genotype 
goi <- c("Il5","Gata2","Cebpb",
         "Retnla",
         "Ear1","Prg2","Prg3", 
         "Siglecf",
         "Il1rl1","Il1r2", "Il10rb","Il1b",
         "Fcgr3",
         "Fosb",
         "Alox5","Alox5ap",
         "Cd47",
         "Spi1",
         "Dennd5a","Dennd1b",
         "Tgfbr1",
         "Rel","Nfkb1")

heatmap_goi_coi(obj, "tissue_cond",goi,"all",length(goi),c(all="#0C6657"),F,T)
