# =======================================================================================
# Purpose: Integration of scRNA-seq and spatial proteomics of the E9 chick cerebellum
# Author: James Li
# Date: July 30, 2024
# 

library(Seurat);library(tidyverse)  
library(SeuratWrappers);library(harmony)

# Load spatial CycIF data and scRNA-seq
seu_sp = readRDS("data/Seurat_aligned_ori.rds")

seu_rna = readRDS("/Volumes/SSD_JamesLi/Experiments1/scRNAseq/Chick_Cb/data/seurat_new.rds")

pdf("figures/UMAP_RNA.pdf", w = 7, h =7)
DimPlot(seu_rna, group.by = "cellType", label = T, repel = T)+NoLegend()+NoAxes()+coord_fixed()
dev.off()

p = FeaturePlot(seu_rna, features = c("FOXP1","CNTNAP5","CRHR2","CNTN5","SULF2","ZFHX4"), combine = F)

for (n in 1:length(p)){
  p[[n]] = p[[n]]+NoLegend()+NoAxes()+coord_fixed()
}

pdf("figures/chick_HCR_probes1.pdf", w = 7, h =7)
cowplot::plot_grid(plotlist = p, ncol = 3)
dev.off()

cell_ids = seu_rna@meta.data %>% 
  filter(cellType %in% c("PC.CRHR2","PC.DPP10","PC.ZFHX4","PC.NR2F2","PC.FOXP4","PC.SPON1")) %>%
  rownames(.)
seu_rna_pc = subset(seu_rna, cells = cell_ids)
table(seu_rna_pc$sample)
table(seu_rna$sample)

Idents(seu_rna_pc) = "cellType"
VlnPlot(seu_rna_pc, features = c("FOXP1","CNTNAP5","CRHR2","CNTN5","SULF2","ZFHX4"),
        pt.size = -1, stack = T,flip = T)

Percentage = table(seu_rna_pc$sample)/table(seu_rna$sample)*100

df = rbind(PC = table(seu_rna_pc$sample), Total = table(seu_rna$sample), Percentage)
write.csv(df, "tables/abundance_PC.csv")


Idents(seu_sp) = "cellClass"
seu_sp_pc = subset(seu_sp, idents = "PC" )
cell_ids = seu_rna@meta.data %>% 
  filter(cellType %in% c("PC.CRHR2","PC.DPP10","PC.FOXP1","PC.NR2F2","PC.FOXP4","PC.SPON1"), sample %in% c("Cb_E7","Cb_E8", "Cb_E9")) %>%
  rownames(.)
seu_rna_pc = subset(seu_rna, cells = cell_ids)
seu_sp_pc$cellType = droplevels(seu_sp_pc$cellType)
seu_rna_pc$cellType = droplevels(seu_rna_pc$cellType)

features_common = intersect(rownames(seu_sp), rownames(seu_rna))
setdiff(rownames(seu_sp), rownames(seu_rna))

anchors_pc <- FindTransferAnchors(reference = seu_rna_pc, query = seu_sp_pc, reduction = 'rpca', normalization.method = "LogNormalize", dims = 1:15, features = features_common)
predicted_pc <- TransferData(anchorset = anchors_pc, refdata = seu_rna_pc$cellType)
seu_sp_pc$RNA_rpca = predicted_pc$predicted.id

CM <- table(seu_sp_pc$RNA_rpca, seu_sp_pc$cellType)
# normalize the matrix with the total cells
CM_norm <- CM/Matrix::rowSums(CM)

pdf("figures/cellType_PConly.pdf", h = 8, w =8.6)
pheatmap::pheatmap(CM_norm, cluster_cols = F, cluster_rows = F,color = ArchR::paletteContinuous("whiteBlue"))
dev.off()

# Transfer the scRNA-seq id to the multplex data
anno = openxlsx::read.xlsx("Annotation_IF1281.xlsx", sheet = "Seurat")
anno_pc = anno %>% filter(cellClass == "PC")
seu_sp_pc$PC_ID = plyr::mapvalues(seu_sp_pc$cellType, anno_pc$cellType, anno_pc$scRNAseq) 

# Show the percentage of different PC subtypes
cell_counts <- seu_sp_pc$PC_ID %>%
  table() %>%                     
  as.data.frame()                 

# Rename the columns 
names(cell_counts) <- c("cellType", "Count")

# Calculate relative abundance
cell_counts$RelativeAbundance <- cell_counts$Count / sum(cell_counts$Count) * 100

# Arrange the cell types by relative abundance in descending order
cell_counts <- cell_counts %>%
  arrange(desc(RelativeAbundance)) %>%
  mutate(cellType = factor(cellType, levels = cellType))  

# Create the bar plot
p1 = ggplot(cell_counts, aes(x = cellType, y = RelativeAbundance, fill = cellType)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative Abundance (%)", title = "Relative Abundance of PC Subtype (CycIF)") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  # Rotate x-axis labels for better visibility

# Print the plot
print(p1)


# Show the percentage of different PC subtypes according to scRNA-seq data
seu_rna_pc$cellType = droplevels(seu_rna_pc$cellType)
cell_count_rna <- seu_rna_pc$cellType %>%
  table() %>%                     
  as.data.frame() 

# Rename the columns 
names(cell_count_rna) <- c("cellType", "Count")
cell_count_rna$cellType = as.character(cell_count_rna$cellType)

# Calculate relative abundance
cell_count_rna$RelativeAbundance <- cell_count_rna$Count / sum(cell_count_rna$Count) * 100

# Arrange the cell types by relative abundance in descending order
cell_count_rna$cellType = factor(cell_count_rna$cellType, levels = levels(cell_counts$cellType))

# Create the bar plot
p2 = ggplot(cell_count_rna, aes(x = cellType, y = RelativeAbundance, fill = cellType)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative Abundance (%)", title = "Relative Abundance of PC Subtype (scRNA-seq)") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  # Rotate x-axis labels for better visibility

# Print the plot
print(p2)

pdf("figures/relative_abundance_new.pdf")
p1
p2
dev.off()


# Show the relative abundance during development according to scRNA-seq
cell_ids1 = seu_rna@meta.data %>% 
  filter(cellType %in% c("PC.CRHR2","PC.DPP10","PC.ZFHX4","PC.NR2F2","PC.FOXP4","PC.SPON1")) %>%
  rownames(.)
seu_rna_pc1 = subset(seu_rna, cells = cell_ids1)
seu_rna_pc1$cellType = as.character(seu_rna_pc1$cellType)
seu_rna_pc1$cellType[seu_rna_pc1$cellType == "PC.ZFHX4"] = "PC.FOXP1"

data <- prop.table(table(seu_rna_pc1$sample, seu_rna_pc1$cellType), margin = 1) %>% as.data.frame() %>% 
  mutate(proportion = Freq*100) %>% 
  rename(stage = Var1, cellType = Var2) %>% 
  mutate(cellType = factor(cellType, levels = levels(cell_counts$cellType)))

# Define a color palette
color_palette <- RColorBrewer::brewer.pal(n = length(unique(data$stage)), name = "Set2")
names(color_palette) <- sort(unique(data$stage))

p3 = ggplot(data, aes(x = cellType, y = proportion, fill = stage)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = color_palette) +
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), legend.position = "top")+
  labs(title = "Relative abundance of PC subtypes per stage", x = "", y = "Proportion (%)")
p3

pdf("figures/RelativeAundancePC.pdf", h =7, w =8)
p3
dev.off()

# Extract the colors
plot_build <- ggplot_build(p1)

# Retrieve colors from the plot build object
used_colors <- plot_build$data[[1]]$fill  # Assuming fill is used in the first layer as it usually is with geom_bar

# To see the unique colors used
unique_colors <- unique(used_colors)
print(unique_colors)


Idents(seu_sp) = "cellType"
clusters <- levels(Idents(seu_sp))
custom_colors <- setNames(rep("lightgrey", length(clusters)), clusters)
custom_colors[as.character(cell_counts$cellType)] <- unique_colors

plot = DimPlot(seu_sp, reduction = "spatial_3d", pt.size = 0.5, split.by = "TR", ncol = 4, 
               group.by = "cellType", cols = custom_colors, label = FALSE, raster = FALSE) +
  NoLegend() + coord_fixed() 

png("figures/PC_only.png", w=2400, h=2500)
print(plot)
dev.off()


data = seu_sp_pc[["spatial_3d"]]@cell.embeddings %>%
  as.data.frame()
data$cellType = seu_sp_pc$PC_ID
data$cellType = as.character(data$cellType)

cluster_slected = c("PC.DPP10",  "PC.CRHR2", "PC.FOXP4",  "PC.NR2F2", "PC.SPON1", "PC.FOXP1")
custom_colors <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33","#FF0303")

data$color = plyr::mapvalues(data$cellType, cluster_slected, custom_colors)

saveRDS(data, "data/for3D.rds")

plot3d( 
  x=data$centroid3D_1, data$centroid3D_2, z=data$centroid3D_3, 
  col = data$color, 
  type = 'p', 
  radius = 0.1,
  aspect = "iso",
  alpha=0.8,
  decorate = FALSE)

rglwidget()

png("figures/3D_Cb_E9.png")

r3dDefaults$windowRect <- c(100, 100, 1000, 1000)
# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 1, 0), rpm = 10), duration = 10 )

# Save like gif
movie3d(
  movie="3dAnimatedScatterplot", 
  spin3d( axis = c(0, 1, 0), rpm = 5),
  duration = 10, 
  dir = "figures/",
  type = "gif", 
  clean = TRUE
)

r3dDefaults$windowRect <- c(100, 100, 1000, 1000)
# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 1, 0), rpm = 10), duration = 10 )

# Save like gif
movie3d(
  movie="3dAnimatedScatterplot", 
  spin3d( axis = c(0, 1, 0), rpm = 5),
  duration = 10, 
  dir = "figures/",
  type = "gif", 
  clean = TRUE
)

library(magick)

# Parameters for the legend
legend_width <- 200
legend_height <- 200
font_size <- 20
line_height <- 30
start_y <- 20

# Create a blank image for the legend
legend_img <- image_blank(width = legend_width, height = legend_height, color = "white")

# Annotate each line with the correct color
for (i in seq_along(cluster_slected)) {
  legend_img <- image_annotate(
    image = legend_img,
    text = cluster_slected[i],            # The cluster name
    color = custom_colors[i],            # Corresponding color
    size = font_size,
    gravity = "northwest",
    location = paste0("+20+", start_y + (i - 1) * line_height)  # Adjust text positioning
  )
}

# Load the GIF
gif <- image_read("figures/3dAnimatedScatterplot.gif")

# Composite the legend onto the GIF
output_gif <- image_composite(gif, legend_img, offset = "+10+10")

# Save the GIF with the legend
image_write(output_gif, "figures/3dAnimatedScatterplot_withLegend.gif")




