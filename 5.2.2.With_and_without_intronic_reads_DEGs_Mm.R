########## This code compares DEGs between intonic and exonic to exonic read mapping only  ##########
### Data from GSE182001

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.3.Functions_gene_expression.R"))

##### Load data 
mouse_il5tg <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_eos_annotation.rds"))
mouse_il5tg_exon_only <- readRDS(file.path(seurat_objects_dir,"Mm_il5tg_eos_forced_cell_determination_exons_only_annotation.rds"))

### Subset cells that are present in both mapping strategies 
intron_and_exon <- rownames(mouse_il5tg@meta.data)
exon_only <- rownames(mouse_il5tg_exon_only@meta.data)

shared_bc <- intersect(intron_and_exon, exon_only)

mouse_il5tg <- subset(mouse_il5tg, cells = shared_bc)
mouse_il5tg_exon_only <- subset(mouse_il5tg_exon_only, cells = shared_bc)

obj <- merge(mouse_il5tg,mouse_il5tg_exon_only)
table(obj$reads)
obj <- JoinLayers(obj)

##### Run DEG analysis and save output csv files 
cell_types <- c("basal","active","circulating","progenitor","immature")
for (i in cell_types) {
  Idents(obj) <- "annotation"
  sub <- subset(obj, idents = i)
  Idents(sub) <- "reads"
  DEG_to_csv_two_cond(sub,"RNA", "intronic_and_exonic","exonic_only",FALSE,0.25,paste0(file.path(gene_mapping_tables_dir, "DEGs_eos_subtype_il5tg_"),i,"_intron_vs_exon_only.csv"))
}

##### Plot the average expression of genes of interest based on the data files produced before 
goi <- c("Myo1d","Elmo1","Asb2","Arhgap15","Atxn1","Syne1","Dock8","Ssh2", # Cytoskeleton 
         "Pla2g4a","Cers6", "Galnt2l",# Lipid metabolism
         "Ikzf3", "Runx1","Morrbid","Ikzf2","Ikzf1","Dach1",# TF
         "Rab8b","Rab2a","Arf4","Dennd1b","Dennd4a","Dennd5a", # Vesicle transport and endocytosis
         "Etv6","Braf","Notch2", "Fli1","Mir142hg", # Hematopoiesis regulation 
         "Cd80","Nfkb1" # Active eos associated
)
obj$anno_cond <- paste0(obj$annotation, "_",obj$reads)
heatmap_goi_coi(obj, "anno_cond",goi,c("Cytoskeleton", "Lipid_metabolism","TF","Vesicle_transport_endocytosis","Heamtopoiesis","Active_eos_associated"), 
                c(8,3,6,6,5,2), c( Cytoskeleton="#0C6657",Lipid_metabolism="#270A7F", TF= "#E293CF",Vesicle_transport_endocytosis="#15E5AE",Heamtopoiesis="#E22B17",
                   Active_eos_associated="#BD7FEA"),F,T)

