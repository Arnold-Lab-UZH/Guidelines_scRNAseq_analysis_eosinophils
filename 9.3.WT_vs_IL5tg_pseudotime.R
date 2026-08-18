########### This code analyses pseudotime trajectory from BM to colon eosinophils comparing il5-tg and WT ##########
### Datasets used: GSE182001, GSE282765

##### Link to libraries and functions
source("0.config.R")
source(file.path(base_dir, "1.1.Packages.R"))

##### read in seurat objects 
obj <- readRDS(file.path(seurat_objects_dir,"wt_il5tg_bm_colon_annotated.rds"))
obj <- JoinLayers(obj)

##### Slingshot analysis 
sce <- sceasy::convertFormat(obj, from="seurat", to="sce")

# plot UMAP grouped by age 
df <- bind_cols(
  as.data.frame(reducedDims(sce)$UMAP.MNN),
  sce$genotype
) %>%
  sample_frac(1)
p1 <- ggplot(df, aes(x = umapmnn_1, y = umapmnn_2, col = ...3)) +
  geom_point(size = .7) +
  scale_color_brewer(palette = "Dark2") +
  labs(col = "genotype")
p1

### Trajectory inference 
sling <- slingshot(sce, clusterLabels = colData(sce)$annotation, reducedDim = "UMAP.MNN"
                   ,omega = TRUE, approx_points = 100)

# Plot trajectory on UMAP 
df <- bind_cols(
  as.data.frame(reducedDims(sling)$UMAP.MNN),
  sling$slingPseudotime_1,
  sling$annotation,
  sling$genotype
) %>%
  sample_frac(1)
colnames(df) <- c("umapmnn_1","umapmnn_2","slingPseudotime_1","annotation","genotype")
curve <- slingCurves(sling)[[1]]
p <- ggplot(df, aes(x = umapmnn_1, y = umapmnn_2, col = slingPseudotime_1)) +
  geom_point(size = .7) +
  scale_color_viridis_c() +
  labs(col = "Pseudotime") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.5) +  theme_classic(base_size = 25) 
ggsave(paste0(file.path(wt_il5tg_plots_dir, "umap_pseudotime.svg")), width = 8, height = 5, plot = p)

### Differential progression between wt and il5tg 
## density plot 
p <- ggplot(df, aes(x = slingPseudotime_1)) +
  geom_density(alpha = .8, aes(fill = genotype), col = "transparent") +
  geom_density(aes(col = genotype), fill = "transparent",
               guide = FALSE, size = 1.5) +
  labs(x = "Pseudotime", fill = "genotype") +
  guides(col = "none", fill = guide_legend(
    override.aes = list(size = 1.5, col = c("#7F7F7C","#8A181A"))
  )) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent")
ggsave(paste0(file.path(wt_il5tg_plots_dir, "pseudotime_density_plot.svg")), width = 8, height = 5, plot = p)

