library(Seurat);library(ggplot2);
library(dplyr)
library(tidyverse)
library(SeuratWrappers)
library(uwot)
adPC <- readRDS("adPC.RDS")

DefaultAssay(adPC) <- "RNA"
ref <- adPC

ref$subcluster = gsub("Purkinje_","", ref$subcluster)
Idents(ref) <- ref$subcluster

ref <- RenameIdents(ref, "Aldoc_1"="Aldoc","Aldoc_2"="Aldoc","Aldoc_3"="Aldoc","Aldoc_4"="Aldoc","Aldoc_5"="Aldoc","Aldoc_6"="Aldoc","Aldoc_7"="Aldoc", "Anti_Aldoc_1"="Anti_Aldoc","Anti_Aldoc_2"="Anti_Aldoc")
DimPlot(ref, label = T)+NoLegend()
prop.table(table(Idents(ref)))
# Aldoc Anti_Aldoc 
# 0.6268486  0.3731514 

ref$Aldoc <- Idents(ref) 

# adPC$cluster = Idents(ref)
# saveRDS(adPC, "adPC.RDS")

ref$Aldoc1 <- ref$subcluster 

ref <- RunUMAP(ref,return.model = T,  umap.method = 'uwot',  reduction = "harmony", dims = 1:30,verbose = FALSE, n.neighbors = 30L,min.dist = 0.3)
PCin <- readRDS("PCin_new.rds")
DimPlot(PCin, group.by = "PC0.8",label = T)+NoLegend()
prop.table(table(PCin$PC0.8))
# PC1        PC2        PC3        PC4        PC5        PC6        PC7        PC8        PC9       PC10       PC11 
# 0.16926202 0.15436696 0.10561950 0.09681787 0.09343263 0.09140149 0.07989167 0.07515234 0.06973595 0.03926879 0.02505078 
# PC2+PC3+PC10 = 0.2992552


DefaultAssay(PCin) <- "RNA"
query <- PCin
anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  query.assay = "RNA",
  normalization.method ="LogNormalize",
  reference.reduction = "pca", k.anchor = 500,
  dims = 1:30,
  recompute.residuals = TRUE
)
saveRDS(anchors, "PC_anchors.rds")

test <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = ref,
  refdata = "Aldoc",
  reference.reduction = "pca", 
  reduction.model = "umap"
)
prop.table(table(test$Aldoc))
# Aldoc Anti_Aldoc 
# 0.91198375 0.08801625 


PCin$Aldoc <- test$predicted.id
DimPlot(PCin, group.by = "Aldoc")

test <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = ref,
  refdata = "Aldoc1",
  reference.reduction = "pca", 
  reduction.model = "umap"
)

PCin$Aldoc1 <- test$predicted.id
DimPlot(PCin, group.by = "Aldoc1")
prop.table(table(test$Aldoc1))
# aAldoc1     aAldoc2      Aldoc2      Aldoc3      Aldoc4      Aldoc5      Aldoc6      Aldoc7 
# 0.008124577 0.094109682 0.023696682 0.029113067 0.098849018 0.437373053 0.244414353 0.064319567 

test <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = ref,
  refdata = "region",
  reference.reduction = "pca", 
  reduction.model = "umap"
)
PCin$region1 <- test$predicted.id
DimPlot(PCin, group.by = "region1")
Idents(PCin) <- "Aldoc1"
levels(PCin) <- c("Aldoc2","Aldoc3","Aldoc4","Aldoc5", "Aldoc6", "Aldoc7","antiAldoc1","antiAldoc2" )
p1 <-  DimPlot(ref, group.by = "Aldoc",label = T)+NoAxes()+NoLegend()
p2 <-  DimPlot(ref, group.by = "region", label = T)+NoAxes()
p3 <-  DimPlot(PCin, group.by = "region1", label = T)+NoAxes()
p4 <-  DimPlot(PCin)+NoAxes()
p5 <-  DimPlot(PCin, group.by = "Aldoc")+NoAxes()

(p1+p2)/(p3+p4+p5)
pdf("adPC_UMAP_alignment3.pdf",w=15,h=10)
(p1+p2)/(p3+p4+p5)
dev.off()

saveRDS(test, "PC_match.rds")
prop.table(table(PCin$))



