library(Seurat);library(tidyverse)
# merging all samples
cdS1 <- readRDS("Multiplexed/integration_P1cKO/data/integrate_celltype.rds")
cdS2 <- readRDS("Multiplexed/Integration_P2-deficient/data/P2-deficient_filtered.rds")
cdS <- merge(cdS1, cdS2)
#feature.common = Reduce(intersect, list(rownames(cdS1),rownames(cdS2)));length(feature.common)
# integrate data
cdS[["RNA"]] <- JoinLayers(cdS[["RNA"]])
cdS[["RNA"]] <- split(cdS[["RNA"]], f = cdS$sample)
cdS <- NormalizeData(cdS) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  # RunPCA(features = feature.common)
  RunPCA(features = rownames(cdS))
cdS <- IntegrateLayers(object = cdS, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = 'harmony', verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:25, min.dist = 0.2) %>%
  FindNeighbors(reduction = "harmony", dims = 1:25, verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE)
p3=DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()
p4=DimPlot(cdS, label = F, group.by = "sample", pt.size = 1)+NoAxes()
p5=DimPlot(cdS, label = F, group.by = "genotype", pt.size = 1)+NoAxes()
p6 = p3+p4+p5
p7=DimPlot(cdS, label = T, split.by = "sample", pt.size = 1)+NoLegend()+NoAxes()
p6/p7
embeddings <- cdS@meta.data[,c("X_centroid","Y_centroid")]
colnames(embeddings) <- c("centroid_1","centroid_2")
cdS[["spatial"]] <- CreateDimReducObject(as.matrix(embeddings), assay = DefaultAssay(cdS), key = "centroid_")
# remove X49531 sample
Idents(cdS) <- "sample"
cdS <- subset(cdS, idents = "X49531", invert=T)
cdS <- NormalizeData(cdS) %>%
FindVariableFeatures() %>%
ScaleData() %>%
# RunPCA(features = feature.common)
RunPCA(features = rownames(cdS))
cdS <- IntegrateLayers(object = cdS, method = HarmonyIntegration, orig.reduction = "pca",
new.reduction = 'harmony', verbose = FALSE) %>%
RunUMAP(reduction = "harmony", dims = 1:25, min.dist = 0.2) %>%
FindNeighbors(reduction = "harmony", dims = 1:25, verbose = FALSE) %>%
FindClusters(resolution = 0.8, verbose = FALSE)
DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()
cdS$res0.8 <- Idents(cdS)
cdS <- FindClusters(cdS,resolution = 1, verbose = FALSE)
cdS$res1 <- Idents(cdS)
Idents(cdS) <- "res1"
p3=DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()
p4=DimPlot(cdS, group.by = "sample", pt.size = 1)+NoAxes()
p5=DimPlot(cdS, group.by = "genotype", pt.size = 1)+NoAxes()
p6 = p3+p4+p5
p7=DimPlot(cdS, label = T, split.by = "sample", pt.size = 1)+NoLegend()+NoAxes()

pdf("figures/UMAP_all.pdf", w =15, h=10)
p6/p7
dev.off()
#Inspect markers and annotate clusters
exp <- AverageExpression(cdS, return.seurat = T)
DoHeatmap(exp, rownames(cdS),label = T,
          draw.lines = F, raster = F, size =2.7, angle = 0)+
  scale_fill_gradientn(colors = colorRampPalette(rev(c("orangered4", "orangered", "white","dodgerblue", "dodgerblue2")))(n = 100))

cdS <- RenameIdents(cdS, c("0"="GC1","1"="IN","2"="PC2","3"="GC2","4"="CN.Gaba","5"="PC3","6"="PC7","7"="PC1","24"="OPC1","8"="OPC2","9"="PC8","10"="GCP4","11"="lCN","12"="mCN","18"="lCN2","13"="GCP1","14"="PC6","15"="UBC", "16"="NPC1","17"="PC9","19"="GCP2.3", "20"="GABA.Pro","21"="PC4","22"="GC3","23"="PCx","25"="Mid.Isl1","26"="cb.bs_border", "27"="GC4","28"="Mid.KO","29"="imPC?","30"="to.filter","31"="PC11","32"="imPCx","33"="bs","34"="bs2","35"="MidO","36"="Mid.X","37"="to.filter","38"="to.filter","39"="to.filter","40"="to.filter","41"="to.filter","42"="to.filter"))

DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()
cdS$celltype_new <- Idents(cdS)
Idents(cdS) <- "celltype_new"
saveRDS(cdS, "cdS_No31.rds")
