########## This code compares DEGs between intonic and exonic to exonic read mapping only  ##########
### Data from GSE282765, Hs PB healthy 

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))
source(file.path(base_dir, "1.3.Functions_gene_expression.R"))

##### Load data 
human_blood <- readRDS(file.path(seurat_objects_dir,"Hs_PB_forced_cell_determination_with_intronic_reads_annotated.rds"))
human_blood_exon_only <- readRDS(file.path(seurat_objects_dir,"Hs_PB_forced_cell_determination_exons_only_annotation.rds"))

### Subset cells that are present in both mapping strategies 
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

##### Run DEG analysis and save output csv files 
conditions <- (as.data.frame(table(obj$phenotype)))$Var1
for (i in conditions) {
  Idents(obj) <- "phenotype"
  sub <- subset(obj, idents = i)
  pb <- AggregateExpression(sub, assays = "RNA", return.seurat = T, group.by = c("condition","reads"))
  Idents(pb) <- "reads"
  DEG_two_cond_pb_DESeq2(pb,"intronic-and-exonic","exonic-only",paste0(file.path(gene_mapping_tables_dir, "DEGs_Hs_blood_"),i,"_intron_vs_exon_only.csv"))
}

##### Plot the average expression of genes of interest based on the data files produced before 
goi <- c("ARHGAP15","ARHGEF18","SSH2","ELMO1","DOCK2","ATXN1","DOCK8", # Actin cytoskeleton
         "SMYD3","PLA2G4A","ACSM3","FAR2","ALOX5", # Lipid metabolism 
         "DACH1","TCF12","RUNX1" ,"RUNX2","IKZF1","IKZF2", # TF
         "DENND1A","DENND4A","RAB2A", # Vesicle Transport and endocytosis 
         "CAMK1D","MARCHF3" # Eosinophil and granulocyte funciton
)

Idents(obj) <- "reads"
p <- DotPlot(obj, features = goi,dot.scale = 10, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave(file.path(gene_mapping_plots_dir, "Blood_PB.svg"), width = 12, height = 6, plot = p)


