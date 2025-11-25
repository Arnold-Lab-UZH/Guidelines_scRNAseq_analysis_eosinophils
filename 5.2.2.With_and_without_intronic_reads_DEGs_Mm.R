########## This code compares DEGs between intonic and exonic to exonic read mapping only  ##########
### Data from GSE182001

##### Set up environment 
setwd("/home/khandl")

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.3.Functions_gene_expression.R")

##### Load data 
mouse_il5tg <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_annotation.rds")
mouse_il5tg_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Mm_il5tg_eos_forced_cell_determination_exons_only_annotation.rds")

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
  DEG_to_csv_two_cond(sub,"RNA", "intronic_and_exonic","exonic_only",FALSE,0.25,paste0("/scratch/khandl/technical/figures/intron_exon/DEGs_eos_subtype_il5tg",i,"_intron_vs_exon_only.csv"))
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
                c(8,3,6,6,5,2),c("#0C6657","#270A7F",  "#E293CF","#15E5AE","#E22B17","#BD7FEA"),
                c( Cytoskeleton="#0C6657",Lipid_metabolism="#270A7F", TF= "#E293CF",Vesicle_transport_endocytosis="#15E5AE",Heamtopoiesis="#E22B17",
                   Active_eos_associated="#BD7FEA"),F,T)

##### DEG analysis 
human_blood <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds")
human_blood_exon_only <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_PB_forced_cell_determination_exons_only_annotation.rds")

### only consider BC that are present in intronic_and_exonic AND exonic_only 
intron_and_exon <- rownames(human_blood@meta.data)
exon_only <- rownames(human_blood_exon_only@meta.data)

shared_bc <- intersect(intron_and_exon, exon_only)

human_blood <- subset(human_blood, cells = shared_bc)
human_blood_exon_only <- subset(human_blood_exon_only, cells = shared_bc)

Idents(human_blood) <- "annotation"
human_blood <- subset(human_blood, idents = "Eosinophils")
Idents(human_blood_exon_only) <- "annotation"
human_blood_exon_only <- subset(human_blood_exon_only, idents = "Eosinophils")

obj <- merge(human_blood,human_blood_exon_only)

obj <- Hs_blood
conditions <- (as.data.frame(table(obj$phenotype)))$Var1
for (i in conditions) {
  Idents(obj) <- "phenotype"
  sub <- subset(obj, idents = i)
  pb <- AggregateExpression(sub, assays = "RNA", return.seurat = T, group.by = c("condition","reads"))
  Idents(pb) <- "reads"
  DEG_two_cond_pb_DESeq2(pb,"intronic-and-exonic","exonic-only",paste0("/scratch/khandl/technical/figures/intron_exon/DEGs_Hs_blood_",i,"_intron_vs_exon_only.csv"))
}

df <- read.csv("/scratch/khandl/technical/figures/intron_exon/DEGs_Hs_blood_healthy_intron_vs_exon_only.csv")

### plot average expression of GOI 
goi <- c("ARHGAP15","ARHGEF18","SSH2","ELMO1","DOCK2","ATXN1","DOCK8", #actin cytoskeleton
         "SMYD3","PLA2G4A","ACSM3","FAR2","ALOX5", #lipid metabolism 
         "DACH1","TCF12","RUNX1" ,"RUNX2","IKZF1","IKZF2", #TF
         "DENND1A","DENND4A","RAB2A", #VESICLE TRANSCPORT and endocytosis 
         "CAMK1D","MARCHF3" #eosinophil and granulocyte funciton
         
)

Idents(obj) <- "reads"
p <- DotPlot(obj, features = goi,dot.scale = 10, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/technical/figures/intron_exon/Blood_PB.svg", width = 12, height = 6, plot = p)


