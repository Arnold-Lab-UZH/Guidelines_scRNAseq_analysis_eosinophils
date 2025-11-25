########## This code compares the quality of additional eosinophils between different empty droplet identification tools ##########
### Datasets used: GSE282765; Hs CRC NAT and tumor;

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.2.Functions_Seurat_integration.R")

##### As a negative control for the quality control extract cells with larger than 25% of mitochondrial reads and lower then 200 features 
###lower 200 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P1",3,0,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,0, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,0,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,0,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,0,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,0,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,0,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,0,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,0,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,0,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,0,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,0,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,0,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl/", "Technical_count_matrices", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,0,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
lower200 <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
lower200 <- JoinLayers(lower200)

### Add mitochondrial percentage per cell 
lower200$percent.mt <- PercentageFeatureSet(lower200, pattern = "^MT-")

### Extract cells with < 200 features 
lower200 <- subset(lower200, subset = nFeature_RNA < 200)
lower200$annotation <- "features_lower200"

### Larger 25% mitochondrial reads 
## P1 
P1_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P1_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST02_Expression_Data.st"), 
                                  "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P1_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST01_Expression_Data.st"), 
                                    "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

## P2 
P2_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P2_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P2_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST08_Expression_Data.st"), 
                                    "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

## P3
P3_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P3_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST09_Expression_Data.st"), 
                                  "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P3_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST10_Expression_Data.st"), 
                                    "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

## P4
P4_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P4_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST04_Expression_Data.st"), 
                                  "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P4_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST03_Expression_Data.st"), 
                                    "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

## P5
P5_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P5_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                  "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P5_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                    "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

## P6
P6_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P6_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST06_Expression_Data.st"), 
                                  "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P6_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST07_Expression_Data.st"), 
                                    "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

## P7 
P7_tumor <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P7_Hs_tumor_forced_cell_determination_intronic_and_exonic_ST12_Expression_Data.st"), 
                                  "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_Hs_data(file.path("/scratch/khandl", "Technical_count_matrices", "P7_Hs_NAT_forced_cell_determination_intronic_and_exonic_ST11_Expression_Data.st"), 
                                    "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")

## Merge samples
larger25 <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))

### Add mitochondrial percentage per cell 
larger25$percent.mt <- PercentageFeatureSet(larger25, pattern = "^MT-")

### Extract cells with larger than 25% genes mapped to mitochondrial reasd 
larger25 <- subset(larger25, subset = percent.mt >= 25)
larger25$annotation <- "mito_larger25"

### Merge larger25 and lower200
merged <- merge(lower200, y = c(larger25),
               add.cell.ids = c("1","2"))
merged <- JoinLayers(merged)
merged <- NormalizeData(merged)

### Save object 
saveRDS( merged, "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_lower200features_larger25mito.rds")

##### Load all Seurat objects 
obj_forced_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
obj_automatic_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_automatic_cell_determination_with_intronic_reads_annotated.rds")
obj_emptyDrops_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_emptyDrops_determination_with_intronic_reads_annotated.rds")
obj_Malat1_all <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_MALAT1_determination_with_intronic_reads_annotated.rds")
obj_with_low200_high25 <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_lower200features_larger25mito.rds")

## Extract Eosinophils and T cells as a control from BD forced 
Idents(obj_forced_all) <- "annotation"
obj_forced <- subset(obj_forced_all, idents = "Eosinophils")
Idents(obj_automatic_all) <- "annotation"
obj_automatic <- subset(obj_automatic_all, idents = "Eosinophils")
Idents(obj_emptyDrops_all) <- "annotation"
obj_emptyDrops <- subset(obj_emptyDrops_all, idents = "Eosinophils")
Idents(obj_Malat1_all) <- "annotation"
obj_Malat1 <- subset(obj_Malat1_all, idents = "Eosinophils")
Tcells <- subset(obj_forced_all, idents = "T")

obj_with_low200_high25$cell_determination <- obj_with_low200_high25$annotation
Tcells$cell_determination <- "T"

##### Extract additional cell IDs from BD forced, MALAT1 threshold and EmptyDrops that are different to BD automatic calls 
automatic_cell_ids <- rownames(obj_automatic@meta.data)
forced_cell_ids <- rownames(obj_forced@meta.data)
emptyDrops_cell_ids <- rownames(obj_emptyDrops@meta.data)
Malat1_cell_ids <- rownames(obj_Malat1@meta.data)

forced_additional <- forced_cell_ids[!forced_cell_ids %in% intersect(automatic_cell_ids, forced_cell_ids)]
emptyDrops_additional <- emptyDrops_cell_ids[!emptyDrops_cell_ids %in% intersect(automatic_cell_ids, emptyDrops_cell_ids)]
Malat1_additional <- Malat1_cell_ids[!Malat1_cell_ids %in% intersect(automatic_cell_ids, Malat1_cell_ids)]

obj_forced <- subset(obj_forced, cells = forced_additional)
obj_emptyDrops <- subset(obj_emptyDrops, cells = emptyDrops_additional)
obj_Malat1 <- subset(obj_Malat1, cells = Malat1_additional)

### Merge all cells togetehr 
obj <- merge(obj_with_low200_high25, y = c(obj_automatic,obj_forced,obj_emptyDrops,obj_Malat1,Tcells), add.cell.ids = c("a","b","c","d","e","f"))
obj <- JoinLayers(obj)

##### Marker gene expression 
markers <- c("CLC","CCR3","ALOX15","ADGRE1","DACH1", "SYNE1","SYNE2", "ADAM8", "MARCHF3","FFAR2","SSH2","PLAUR","CCL4L2","AREG","NFKB1","RELB","NCF2","CAMK1D", "OLR1","CD3E","ICOS","CD8A","CD4")

Idents(obj) <- "cell_determination"
p <- DotPlot(obj, features = markers,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits =c("features_lower200","mito_larger25","T","MALAT1" ,"EmptyDrops", "forced","automatic"))
ggsave("/scratch/khandl/technical/figures/summary_empty/DotPlot_markers_real_different_approaches.svg", width = 15, height = 6, plot = p)

p <- VlnPlot(obj, features= "nFeature_RNA", group.by = "cell_determination", pt.size = 0, cols = c("#7251A1","#6DA0D5" ,"#6D7272","#F4C914" ,"#23803F", "#020202","#5BC8D9" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  +  
  scale_x_discrete(limits =c("automatic", "forced","EmptyDrops","MALAT1" ,"T","mito_larger25","features_lower200"))
print(p)
ggsave("/scratch/khandl/technical/figures/summary_empty/nFeature.svg", width = 8, height = 8, plot = p)

p <- VlnPlot(obj, features= "percent.mt", group.by = "cell_determination", pt.size = 0, cols = c("#7251A1","#6DA0D5" ,"#6D7272","#F4C914" ,"#23803F", "#020202","#5BC8D9" ) )+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis()  + 
  scale_x_discrete(limits =c("automatic", "forced","EmptyDrops","MALAT1" ,"T","mito_larger25","features_lower200"))
print(p)
ggsave("/scratch/khandl/technical/figures/summary_empty/percent_mt.svg", width = 8, height = 8, plot = p)

