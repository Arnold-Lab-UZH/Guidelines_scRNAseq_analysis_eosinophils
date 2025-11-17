########## This code quantifies ambient RNA content  ##########

##### Set up environment 
setwd("/home/khandl")

##### link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

# GSE282765 and GSE182001
# take forced cell determination and intronic + exonic reads

##### read R objects 
obj <- readRDS( "/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")

current.cluster.ids <- c("?","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                         "Monocytes","Neutrophils","PCs","TAMs", "T")
new.cluster.ids <- c("Undefined","B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "lowQ", "Macrophages","Mast","Mixed",
                     "Monocytes","Neutrophils","PCs","TAMs", "T")
obj$annotation <- plyr::mapvalues(x = obj$annotation, from = current.cluster.ids, to = new.cluster.ids)

##### run for each batch 
### batch 1 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp1")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P1_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.04845,0.04343,0.06611,0.04082,0.07385,0.05874,0.0516,0.04124,0.05628,0.05596,0.27236,0.06338,0.04759)
df1 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P1_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.05632,0.0521,0.08219,0.0454,0.07586,0.0818,0.05453,0.04784,0.0476,0.05614,0.25834,0.05595,0.04838)
df2 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch1_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 2 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp2")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
# NAT 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P2_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(0.03929,0.03781,0.05187,0.04131,0.07206,0.0538,0.03272,0.03077,0.0273,0.057,0.25303,0.0,0.03631)
df3 <- data.frame(cell_types, contamination_values)

# tumor 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P2_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.04007,0.0336,0.05813,0.03647,0.04717,0.03025,0.0285,0.04042,0.03707,0.05955,0.25303,0.02641,0.03442)
df4 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch2_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 3 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp3")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P3_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.06318,0.04467,0.05507,0.05672,0.05602,0.05571,0.04531,0.04723,0.04045,0.05533,0.26362,0.04503,0.05564)
df5 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P3_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.05935,0.04418,0.05518,0.05808,0.0365,0.04387,0.04012,0.04586,0.04517,0.08653,0.29625,0.03249,0.05071)
df6 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch3_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 4
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp4")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P4_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.08528,0.04506,0.0663,0.06667,0.05886,0.05786,0.04858,0.04507,0.04935,0.04377,0.20537,0.02836, 0.05982)
df7 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P4_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.07036,0.05162,0.06212,0.06467,0.04333,0.04975,0.0379,0.0521,0.05027,0.04917,0.19097,0.03197,0.06437)
df8 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch4_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 5
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp5")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P5_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.02404,0.02201,0.01735,0.01237,0.03295,0.01777,0.02067,0.02075,0.02087,0.00928,0.17625,0.01571,0.02276)
df9 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P5_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.01996,0.02087,0.02113,0.01198,0.0401,0.01948,0.01792,0.01949,0.01876,0.00967,0.19891,0.01607,0.01863)
df10 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch5_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 6 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp7")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P6_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.01931, 0.0222,0.01625, 0.0053,0.02868,0.01294,0.01964,0.01846,0.02151,0.01852,0.2251,0.02145,0.01658)
df11 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P6_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.0232,0.02355,0.0209,0.006,0.02362,0.02025,0.02132,0.02115,0.01867,0.03853,0.23419,0.02255,0.01855)
df12 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch6_corrected = ContaminationCorrection(sub,rownames(GCGs))

### batch 8 
Idents(obj) <- "experiment"
sub <- subset(obj, idents = "Exp8")

## define contamination-causing genes (GCGs) 
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "annotation"
GCGs <-  ContaminationDetection(sub,out_path.plot = "/scratch/khandl/technical/figures/Ambient_RNA/")
rownames(GCGs)

## quantify contamination per celltype 
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")

# NAT 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P7_tissue_ctrl")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
cell_types <- c("B","DCs","Endothelial", "Eosinophils","Epithelial","Fibroblasts", "Macrophages","Mast",
                "Monocytes","Neutrophils","PCs","TAMs", "T")
contamination_values <- c(0.05741,0.05261,0.05691,0.07853,0.0,0.05717,0.05546,0.05854,0.05437,0.05564,0.20399,0.0,0.05951)
df13 <- data.frame(cell_types, contamination_values)

# tumor 
Idents(sub) <- "condition"
sub_2 <- subset(sub, idents = "P7_tumor")
for(i in cell_types) {
  Idents(sub_2) <- "annotation"
  sub2 <- subset(sub_2, idents = i)
  ContaminationQuantification(sub2,rownames(GCGs))
  print(i)
}
contamination_values <- c(0.06375,0.04987,0.05036,0.06502,0.03774,0.05851,0.04837,0.0534,0.05087,0.05048,0.21745,0.03911,0.05703)
df14 <- data.frame(cell_types, contamination_values)

## correct the contamination
Idents(sub) <- "annotation"
batch7_corrected = ContaminationCorrection(sub,rownames(GCGs))

