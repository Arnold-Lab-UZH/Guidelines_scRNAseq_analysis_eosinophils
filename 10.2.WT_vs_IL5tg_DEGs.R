########## This code quantifies ambient RNA content using SoupX  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.3.Functions_gene_expression.R")

# GSE282765 
# take forced cell determination and intronic + exonic reads

##### read in seurat objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/wt_il5tg_bm_colon_annotated.rds")
obj <- JoinLayers(obj)

##### per tissue 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "colon")
Idents(sub) <- "genotype"
DEG_to_csv_two_cond(sub,"RNA", "il5tg","wt",FALSE,0.25,"/scratch/khandl/technical/figures/wt_il5tg/DEGs_colon_il5tg_vs_wt.csv")

sub <- subset(obj, idents = "bm")
Idents(sub) <- "genotype"
DEG_to_csv_two_cond(sub,"RNA", "il5tg","wt",FALSE,0.25,"/scratch/khandl/technical/figures/wt_il5tg/DEGs_bm_il5tg_vs_wt.csv")


obj$tissue_cond <- paste0(obj$tissue, "_",obj$condition)

##### plot average expression of GOI between tissue and genotype 
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

heatmap_goi_coi(obj, "tissue_cond",goi,"all",length(goi),"#0C6657",c(all="#0C6657"),F,T)

heatmap_goi_coi(obj, "anno_cond",goi,c("cytoskeleton", "lipid_metabolism","TF","Vesicle_transport_endocytosis","heamtopoiesis","active_eos_associated"), 
                c(8,3,6,6,5,2),c("#0C6657","#270A7F",  "#E293CF","#15E5AE","#E22B17","#BD7FEA"),
                c( cytoskeleton="#0C6657",lipid_metabolism="#270A7F", TF= "#E293CF",Vesicle_transport_endocytosis="#15E5AE",heamtopoiesis="#E22B17",
                   active_eos_associated="#BD7FEA"),F,T)
