########### This code compares wilcox, MAST and DESeq2 statistical approaches for DEGs in eosinophils against macrophages ##########
### Datasets used: GSE282765

##### Set up environment 
setwd("/home/khandl")

##### Link to libraries and functions
source("~/Projects/Guidelines_scRNAseq_analysis_eosinophils/1.1.Packages.R")

##### Load data 
obj <- readRDS("/scratch/khandl/technical/seurat_objects/Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds")
Idents(obj) <- "reads"
obj <- subset(obj, idents = "intronic_and_exonic")
Idents(obj) <- "cell_determination"
obj <- subset(obj, idents = "forced")

### wilcox
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindMarkers(object = obj, ident.1 = "Eosinophils",ident.2 ="Macrophages", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_wilcox_Eos_vs_Mac.csv")

### MAST
markers <- FindMarkers(object = obj,  ident.1 = "Eosinophils",ident.2 ="Macrophages", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data",test.use = "MAST", latent.vars = c("nFeature_RNA","nCount_RNA"))
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_MAST_Eos_vs_Mac.csv")

### DESeq2 
pb <- AggregateExpression(obj, assays = "RNA", return.seurat = T, group.by = c("condition","annotation"))
Idents(pb) <- "annotation"
markers <- FindMarkers(object = pb, ident.1 = "Eosinophils",ident.2 = "Macrophages", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data",test.use = c("DESeq2"))
write.csv(markers, "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_DESeq2_Eos_vs_Mac.csv")

### Compare Eosinophil markers 
wilcox <- read.csv( "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_wilcox_Eos_vs_Mac.csv")
MAST <- read.csv( "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_MAST_Eos_vs_Mac.csv")
DESeq2 <- read.csv( "/scratch/khandl/technical/figures/Eos_marker/Eos_marker_Hs_NAT_tumor_GSE282765_DESeq2_Eos_vs_Mac.csv")

wilcox <- wilcox[wilcox$p_val_adj <= 0.05 & wilcox$avg_log2FC >0,]
wilcox_top10 <- wilcox[order(-wilcox$avg_log2FC), ][1:10, ]
MAST <- MAST[MAST$p_val_adj <= 0.05 & MAST$avg_log2FC >0,]
MAST_top10 <- MAST[order(-MAST$avg_log2FC), ][1:10, ]
DESeq2 <- DESeq2[DESeq2$p_val_adj <= 0.05 & DESeq2$avg_log2FC >0,]
DESeq2_top10 <- DESeq2[order(-DESeq2$avg_log2FC), ][1:10, ]

x <- list("wilcox" = wilcox$X, "MAST" =MAST$X,"DESeq2"= DESeq2$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#EABBC4")

x <- list("wilcox" = wilcox_top10$X, "MAST" =MAST_top10$X, "DESeq2"=DESeq2_top10$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) + scale_fill_gradient(low = "white", high = "#EABBC4")

shared <- Reduce(intersect, list(wilcox_top10$X,MAST_top10$X,DESeq2_top10$X))

# merge tables 
wilcox_top10 <- wilcox_top10[, colnames(wilcox_top10) %in% c("X","avg_log2FC","p_val_adj")]
wilcox_top10$method <- "Wilcox"
MAST_top10 <- MAST_top10[, colnames(MAST_top10) %in% c("X","avg_log2FC","p_val_adj")]
MAST_top10$method <- "MAST"
DESeq2_top10 <- DESeq2_top10[, colnames(DESeq2_top10) %in% c("X","avg_log2FC","p_val_adj")]
DESeq2_top10$method <- "DESeq2"

df <- rbind(wilcox_top10,MAST_top10)
df <- rbind(df,DESeq2_top10)

### visualize p-value zero inflation 
genes_to_label <- shared
p <- ggplot(DESeq2, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = X %in% genes_to_label)) + 
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(
    data = subset(DESeq2, X %in% genes_to_label),
    aes(label = X),
    size = 4
  ) + ylim(0,300) + 
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
ggsave("/scratch/khandl/technical/figures/Eos_marker/DESeq2_scatter.svg", width = 10, height = 8, plot = p)

p <- ggplot(wilcox, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = X %in% genes_to_label)) + 
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(
    data = subset(DESeq2, X %in% genes_to_label),
    aes(label = X),
    size = 4
  ) + ylim(0,300) + 
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
ggsave("/scratch/khandl/technical/figures/Eos_marker/wilcox_scatter.svg", width = 10, height = 8, plot = p)

p <- ggplot(MAST, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = X %in% genes_to_label)) + 
  scale_color_manual(values = c("black", "red")) +
  geom_text_repel(
    data = subset(DESeq2, X %in% genes_to_label),
    aes(label = X),
    size = 4
  ) + ylim(0,300) + 
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
ggsave("/scratch/khandl/technical/figures/Eos_marker/MAST_scatter.svg", width = 10, height = 8, plot = p)





