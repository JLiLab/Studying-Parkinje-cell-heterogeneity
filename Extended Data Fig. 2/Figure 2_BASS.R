# =======================================================================================
# Purpose: Cell clusering and domain segmentation based on imputed spatial transcriptome
# Author: James Li
# Date: Jan 05, 2025
# 

library(Seurat)
library(BASS)
library(tidyverse)

obj = readRDS("data/merge_4x.rds")

obj.list = SplitObject(obj, split.by = "ID")

# Initialize lists for gene expression and spatial coordinates
gene_expression_list <- list()
coordinates_list <- list()

id_selected = setdiff(names(obj.list), c("53168_TR11", "53168_TR12"))

# Loop through each section
for (id in id_selected) {
  # Store results in the lists
  gene_expression_list[[id]] <- GetAssayData(obj.list[[id]], layer = "counts", assay = "RNA")
  coordinates_list[[id]] <- as.matrix(obj.list[[id]]@meta.data[, c("centroid_1", "centroid_2")])
}

rm(obj)

BASS_obj <- createBASSObject(
  X = gene_expression_list, 
  xy = coordinates_list, 
  C = 40, 
  R = 30, 
  beta_method = "SW"
)

BASS_obj <- BASS.preprocess(BASS_obj)
BASS_obj <- BASS.run(BASS_obj)
BASS_obj <- BASS.postprocess(BASS_obj)

png("figures/BASS_tracePlot.png", w = 500, h =300)
plot(1:BASS_obj@burnin, BASS_obj@samples$beta,xlab = "Iteration", 
     ylab = expression(beta), type = "l")
dev.off()

obj$cluster_BASS1 <- factor(unlist(BASS_obj@results$c))
obj$domain_BASS1 <- factor(unlist(BASS_obj@results$z))

dir.create("figures/cluster_inspection_BASS1")
Idents(obj) = "cluster_BASS1"
for (n in levels(Idents(obj))){
  cells = WhichCells(obj, idents = n)
  p = DimPlot(obj, reduction = "spatial", split.by = "ID", cells.highlight = cells, ncol = 4, raster = F, pt.size = 0.5)+
    NoLegend()+NoAxes()+coord_fixed()+ggtitle(paste0("cluster_",n))
  png(paste0("figures/cluster_inspection_BASS1/cluster_",n,".png"), w = 3700, h =800)
  print(p)
  dev.off()
}

dir.create("figures/domain_BASS1")
Idents(obj) = "domain_BASS1"
for (n in levels(Idents(obj))){
  cells = WhichCells(obj, idents = n)
  p = DimPlot(obj, reduction = "spatial", split.by = "ID", cells.highlight = cells, ncol = 4, raster = F, pt.size = 0.5)+
    NoLegend()+NoAxes()+coord_fixed()+ggtitle(paste0("cluster_",n))
  png(paste0("figures/domain_BASS1/cluster_",n,".png"), w = 3700, h =800)
  print(p)
  dev.off()
}

saveRDS(obj, "data/merge_4x.rds")

# Detect marker genes for each domains
library(COSG)
Idents(seu) ="domain_BASS1"
marker_cosg <- cosg(seu, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=50)

# Extract gene markers, scores, and clusters
gene_markers <- marker_cosg$names  # List of genes for each cluster
gene_scores <- marker_cosg$scores  # Corresponding scores for each gene
clusters <- names(gene_markers)    # Cluster names (or group labels)

# Create a data frame by combining genes, scores, and clusters
df_cosg_markers <- data.frame(
  gene = unlist(gene_markers),                  # Unlist all gene names into one column
  score = unlist(gene_scores),                  # Unlist all scores into one column
  cluster = rep(clusters, sapply(gene_markers, length))  # Repeat cluster names based on the number of genes per cluster
)

df_cosg_markers$TF = mGenes$Type[match(df_cosg_markers$gene, mGenes$mgi_symbol)]
df_cosg_markers$description = mGenes$description[match(df_cosg_markers$gene, mGenes$mgi_symbol)]

head(df_cosg_markers)
openxlsx::write.xlsx(df_cosg_markers,asTable = TRUE, "tables/COSG_Markers_domain_BASS1.xlsx")

# Annotate domains manually in excel and read in the information 
anno = openxlsx::read.xlsx("tables/COSG_Markers_domain_BASS1.xlsx")
anno = anno[!duplicated(anno$cluster),]
seu$cellGroup = plyr::mapvalues(seu$domain_BASS1, anno$cluster, anno$cellType)

Idents(seu) = "cellGroup"
clusters <- levels(Idents(seu))

# Initialize all clusters with a base color (e.g., gray90)
custom_colors <- setNames(rep("gray90", length(clusters)), clusters)
cluster_selected <- c("GCP", "CN.gaba", "FN", "IP", "DN", "PC1", "PC2", "PC7", "PC6", "PC9","GABA.Pre","UBC")
custom_colors[cluster_selected] <- RColorBrewer::brewer.pal(n = length(cluster_selected), name = "Paired")

# Select four representative sections for the Figure
Idents(seu) = "ID"
obj = subset(seu, idents = c("43430_TR8","48793_TR7","53168_TR11","56831_TR29"))

pdf("figures/Spatial_domain3.pdf", w = 12, h =6)
DimPlot(obj, reduction = "spatial", pt.size = 0.1, group.by = "cellGroup",split.by = "ID",raster = F, ncol=2,
        cols = custom_colors)+NoAxes()+NoLegend()+coord_fixed()+ggtitle("")
dev.off()


# Print color pallets for each domains
filtered_colors <- custom_colors[custom_colors != "gray90"]

pdf("figures/Spatial_domain_pellet.pdf", w = 5, h =3)
barplot(rep(1, length(filtered_colors)), col = filtered_colors, names.arg = names(filtered_colors), las = 2)
dev.off()


# Examine the consistency of cell abundance BASS-defined domains
seu_sub = subset(seu, idents = cluster_selected)
seu_sub$cellGroup = droplevels(seu_sub$cellGroup)
df = prop.table(table(seu_sub$cellGroup, seu_sub$sample), margin = 2) %>% 
  as.data.frame() %>% 
  set_names("cellGroup","sample","ratio") %>% 
  mutate(percentage = ratio * 100)

df$cellGroup = droplevels(df$cellGroup)
df$sample = as.character(df$sample)

df$cellGroup = as.character(df$cellGroup)
df$cellGroup = factor(df$cellGroup, levels = gtools::mixedsort(unique(df$cellGroup)))
df$cellGroup = factor(df$cellGroup, 
                      levels = c("GCP", "FN", "IP", "DN", "UBC","CN.gaba","PC1", "PC2", "PC7", "PC6", "PC9","GABA.Pre","NPC"))

pdf("figures/composition_BASS_cellGroup.pdf", h = 5, w = 5)
ggplot(df, aes(x = sample, y = percentage, fill = sample))+
  geom_bar(stat = "identity")+
  facet_wrap(~cellGroup,ncol = 4)+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_blank(), legend.position = "top") +
  ylab("Percentage (%)")+xlab("")
dev.off()  


