########### This code analyses macrophages and eosinophils in terms of their marker gene expression - single- vs. multi- genes ##########
### Datasets used: GSE282765

##### link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### load data
obj <- readRDS(file.path(seurat_objects_dir,"Hs_tumor_NAT_forced_cell_determination_with_intronic_reads_annotated.rds"))
Idents(obj) <- "annotation"
sub <- subset(obj, idents = c("Eosinophils","Macrophages"))
#Identify the top gene expressed in both clusters 
Idents(sub) <- "annotation"
markers <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", layer = "data", test.use = "MAST")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
top10 <- (as.data.frame(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)))$gene

### DotPlot of top DEGs 
p <- DotPlot(sub, features = top10,dot.scale = 15, scale = FALSE, assay = "RNA",cols = c("white","darkred")) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 
ggsave(file.path(dropout_plots_dir, "Mac_Eos_marker.svg"), width = 15, height = 6, plot = p)

##### plot the percentage of single- and multiple marker genes 
DefaultAssay(sub) <- "RNA"
sce <- as.SingleCellExperiment(sub)

markers <- list(Mac_C1QA = c("C1QA"),
                Mac_C1QC = c("C1QC"),
                Eos = c("FFAR2","MARCHF3","DACH1","SYNE1","CLC"),
                Eos_FFAR2 = c("FFAR2"),
                Eos_SYNE1 = c("SYNE1"),
                Eos_MARCHF3 = c("MARCHF3"),
                Eos_DACH1 = c("DACH1"),
                Eos_CLC = c("CLC")
                )

colData(sce)$cluster <- sce$annotation
df <- plotDecontXMarkerPercentage_df(sce, markers = markers,assayName = "counts", z = "annotation")

## Mac markers 
df1 <- df[df$markerLabels %in% "Mac_C1QA",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#82C341", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Mac_marker_C1QA.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Mac_C1QC",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#82C341", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Mac_marker_C1QC.svg"), width = 10, height = 8, plot = p)

## Eosinohil markers 
df1 <- df[df$markerLabels %in% "Eos_FFAR2",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_FFAR2.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Eos_DACH1",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_DACH1.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Eos_CLC",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_CLC.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Eos_MARCHF3",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_MARCHF3.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Eos_SYNE1",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_SYNE1.svg"), width = 10, height = 8, plot = p)

df1 <- df[df$markerLabels %in% "Eos",]
p <- ggplot(df1, aes(x = cellType, y =  percent, fill = assay)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  geom_text(aes(label = round(percent, 1)),  # round to 1 decimal
            position = position_dodge(width = 0.8),  # same dodge as bars
            vjust = -0.5,   # slightly above the bar
            size = 3) +     # text size 
  scale_fill_manual(values = c("counts" = "#E23128", name = "Assay")) +  # manual colors & legend title
  theme_bw() + ylim(0,100) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave(file.path(dropout_plots_dir, "Eos_marker_all5.svg"), width = 10, height = 8, plot = p)
