library(Seurat);library(dplyr);library(Matrix);library(SeuratWrappers);library(ggplot2);library(patchwork);library(tidyverse);
library(reshape2);library(scales);
cdS <-readRDS("E16_18_All.rds")
Idents(cdS) <- "data" 
cdS <- subset(cdS, idents = c("In-house"))
cdS[["percent.mito"]] <- PercentageFeatureSet(object = cdS, pattern = "^mt-")
Idents(cdS) <- "dataset"
E16i <- subset(cdS, idents = "E16i")
E18i <- subset(cdS, idents = "E18i")
E16n <- subset(cdS, idents ="E16n")
VlnPlot(object = cdS,group.by = "dataset", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)
VlnPlot(object = E16i,group.by = "dataset", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)
E16i <- subset(E16i, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mito < 8 & nCount_RNA < 20000)
E16i <- NormalizeData(E16i) %>% ScaleData() %>%FindVariableFeatures()

VlnPlot(object = E18i,group.by = "dataset", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)
E18i <- subset(E18i, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mito < 7.5 & nCount_RNA < 20000)
E18i <- NormalizeData(E18i) %>% ScaleData() %>%FindVariableFeatures()

VlnPlot(object = E16n,group.by = "dataset", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0)
E16n <- subset(E16n, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mito < 15 & nCount_RNA < 30000)
E16n <- NormalizeData(E16n) %>% ScaleData() %>%FindVariableFeatures()

cdS <- merge(E16i,E18i)
cdS <- merge(cdS,E16n)
p1 <- FeatureScatter(object = cdS, feature1 = "nCount_RNA", feature2 = "percent.mito",group.by = "dataset")
p2 <- FeatureScatter(object = cdS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "dataset")
cowplot::plot_grid(p1,p2)

# Use FastMNN to integrate data
cdS <- RunFastMNN(object.list = SplitObject(cdS, split.by = "dataset"))
cdS <- RunUMAP(cdS, reduction = "mnn", dims = 1:30, reduction.name = "umap",n.neighbors = 35L,min.dist = 0.4,umap.method = "uwot",return.model = T)
cdS <- FindNeighbors(cdS, reduction = "mnn", dims = 1:30)
cdS <- FindClusters(object = cdS, verbose = FALSE,resolution = 1)
p1 <- DimPlot(cdS, label = T)+NoAxes()+NoLegend()             
p2 <- DimPlot(cdS, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(cdS, group.by = "stage")+NoAxes()
p4 <- DimPlot(cdS, group.by = "dataset")+NoAxes()
(p1+p2)/(p3+p4)
DimPlot(cdS, label = T)+NoAxes()+NoLegend()

pdf("figures/cdS_mnn.pdf", h =10, w = 10)
(p1+p2)/(p3+p4)
dev.off()

pdf("figures/cdS_UMAP.pdf", h =10, w = 10)
p1
dev.off()

cdS <- RunTSNE(cdS,reduction = "mnn",dims = 1:50)
DimPlot(cdS, label = T,reduction = "tsne")+NoAxes()+NoLegend()

saveRDS(cdS, "data/cdS.rds")


