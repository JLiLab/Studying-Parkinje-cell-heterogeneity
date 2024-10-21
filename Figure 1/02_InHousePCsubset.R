library(Seurat);library(dplyr);library(Matrix);library(SeuratWrappers);library(ggplot2);library(patchwork);library(tidyverse);
library(reshape2);library(scales);
cdS <- readRDS("data/cdS.rds")
# subset PC and CN
PC.CN <-  subset(cdS, idents = c("8","18","21"))
PC.CN <- NormalizeData(PC.CN) %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
# Use FastMNN to integrate the data
PC.CN <- RunFastMNN(object.list = SplitObject(PC.CN, split.by = "dataset"))
PC.CN <- RunUMAP(PC.CN, reduction = "mnn", dims = 1:20, reduction.name = "umap",n.neighbors = 20L)
PC.CN <- FindNeighbors(PC.CN, reduction = "mnn", dims = 1:30)
PC.CN <- FindClusters(object = PC.CN, verbose = FALSE, resolution = 0.8)

p1 <- DimPlot(PC.CN, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(PC.CN, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(PC.CN, group.by = "stage")+NoAxes()
p4 <- DimPlot(PC.CN, group.by = "dataset")+NoAxes()
(p1+p2)/(p3+p4)

FeaturePlot(PC.CN, "Pax2",label = T)
pdf("figures/PC.CN_mnn.pdf", h =10, w = 10)
(p1+p2)/(p3+p4)
dev.off()

p1
# subset PC only
FeaturePlot(PC.CN, c("Foxp2","Sox14","Meis2","Tfap2b"))
PC <- subset(PC.CN, idents = c("3","10"),invert=T)    
PC <- NormalizeData(PC) %>% ScaleData()  %>% FindVariableFeatures() %>% RunPCA()
# Use FastMNN to integrate the data
PC <- RunFastMNN(object.list = SplitObject(PC, split.by = "dataset"))
PC <- RunUMAP(PC, reduction = "mnn", dims = 1:30, reduction.name = "umap",n.neighbors = 30L)
PC <- FindNeighbors(PC, reduction = "mnn", dims = 1:30)
PC <- FindClusters(object = PC, verbose = FALSE, resolution = 1)

p1 <- DimPlot(PC, label = T)+NoAxes()+NoLegend()
p2 <- DimPlot(PC, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(PC, group.by = "stage")+NoAxes()
p4 <- DimPlot(PC, group.by = "dataset")+NoAxes()
(p1+p2)/(p3+p4)

pdf("figures/PC_mnn_res0.8.pdf", h =10, w = 10)
(p1+p2)/(p3+p4)
dev.off()


saveRDS(PC.CN, "data/PC.CN.rds")
saveRDS(PC, "data/PC.rds")