########### This code analyses granulogenesis and antimicrobial peptide score  ##########
### Datasets used: GSE182001, GSE282765

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### read in seurat objects 
obj <- readRDS(file.path(seurat_objects_dir,"wt_il5tg_bm_colon_annotated.rds"))
obj <- JoinLayers(obj)

##### BM 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "bm")

### Signature scores 
Granulogenesis_list <- list(c("Prg2","Prg3",  "Epx", "Ear6", "Ear1", "Ear2"))
sub <-AddModuleScore(sub, features= Granulogenesis_list,name = "Granulogenesis")
Antimicrobial_peptides <- list(c("S100a8","S100a9","S100a11","S100a12","Ltf","Camp","Lcn2"))
sub <-AddModuleScore(sub, features= Antimicrobial_peptides,name = "Antimicrobial_peptides")

p <- VlnPlot(sub, features= "Granulogenesis1", group.by = "genotype", pt.size = 1) + stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(file.path(wt_il5tg_plots_dir, "Granulogenesis_bm.svg"), width = 5, height = 8, plot = p)

p <- VlnPlot(sub, features= "Antimicrobial_peptides1", group.by = "genotype", pt.size = 1) + stat_compare_means(method = "wilcox.test", label = "p.format")
ggsave(file.path(wt_il5tg_plots_dir, "Antimicrobial_bm.svg"), width = 5, height = 8, plot = p)
