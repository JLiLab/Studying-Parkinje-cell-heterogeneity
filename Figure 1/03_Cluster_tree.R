
library(Seurat);library(tidyverse);library(clustree)

PCin <- readRDS("../data/PCin_new.rds")

seu = PCin
seu@meta.data = seu@meta.data[,c(1:5,40)]
DefaultAssay(seu) <- "RNA"
for (res in c(0,0.1,0.3,0.5,0.7,0.9,1.1)){
  seu <- FindClusters(seu, resolution = res)
}

df = seu@meta.data

# add cell type information
df$RNA_snn_res.0.8 = df$PC0.8

pdf("../figures/clusterTree_PC.pdf", w =8, h=8)
clustree(df, prefix = "RNA_snn_res.")
dev.off()

old = c(0:10)
new = c("PC1","PC2","PC4","PC3","PC6","PC7","PC8","PC9","PC5","PC10","PC11")
seu$PC= seu$RNA_snn_res.0.7 = plyr::mapvalues(seu$RNA_snn_res.0.7, old, new)

df = seu@meta.data

pdf("../figures/clusterTree_PC1.pdf", w =7, h=7)
clustree(df, prefix = "RNA_snn_res.")
dev.off()


library(COSG)
Idents(seu) = seu$PC
marker_cosg <- cosg(seu, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=3)
all_markers<- marker_cosg$names %>%
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

p <- DotPlot(object = seu, features = all_markers)+ coord_flip()+ RotatedAxis()
df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)

p1 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 

pdf("../figures/PC_markers_newCluster.pdf", h =6, w=5)
p1
dev.off()

saveRDS(seu, "../data/PC_newCluster.rds")
