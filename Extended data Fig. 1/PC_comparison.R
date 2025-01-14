# =======================================================================================
# Goal: Examine the robustness of PC subtype classification using publicly available data
# Author: Nagham Khouri-Farah1
# Date: Sep 1, 2023
# 

library(Seurat);
library(tidyverse);

E1618a <- readRDS("/data/E16_18_All.rds")
DimPlot(E1618a,group.by = "data")
Idents(E1618a) <- "data"

# Extract Taylor's data
E1618t <- subset(E1618a, idents="Taylor")
DimPlot(E1618t, group.by = "cellType1", label = T)+NoLegend()+NoAxes()

E1618t <- NormalizeData(E1618t) %>% ScaleData() %>%FindVariableFeatures() %>%
  RunPCA(verbose = FALSE)
E1618t <- RunUMAP(E1618t, dims = 1:50, reduction.name = "umap",n.neighbors = 35L)
E1618t <- FindNeighbors(E1618t,  dims = 1:50)
E1618t <- FindClusters(object = E1618t, verbose = FALSE)
p1 <- DimPlot(E1618t, label = T)+NoAxes()+NoLegend()             
p2 <- DimPlot(E1618t, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(E1618t, group.by = "stage")+NoAxes()
p4 <- DimPlot(E1618t, group.by = "dataset")+NoAxes()
(p1+p2)/(p3+p4)
DimPlot(E1618t, label = T,)+NoAxes()+NoLegend()

#subset PC, IN and CN 
CNPCt <- subset(E1618t, idents=c("12","15","8"))
CNPCt <- NormalizeData(CNPCt) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)
CNPCt <- RunUMAP(CNPCt,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.3) %>%
  FindNeighbors( dims = 1:25, verbose = FALSE) %>%
  FindClusters(resolution = 1.5, verbose = FALSE)
DimPlot(CNPCt, label = T)+NoAxes()+NoLegend()
FeaturePlot(CNPCt, c("Foxp2","Meis2","Pax2","Pax6"))

#remove non-PC cells
PCt <- subset(CNPCt, idents=c(1,5,8,12),invert=T)
PCt <- NormalizeData(PCt) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE) 
PCt <- RunUMAP(PCt,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.3)%>%
  FindNeighbors( dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 2, verbose = FALSE)
DimPlot(PCt, label = T)+NoAxes()+NoLegend()


# Extract La Manno's data from the integrated dataset
E1618m <- subset(E1618a, idents="LaManno")
EmbeddingsMA <- E1618m@reductions$umap@cell.embeddings
DimPlot(E1618m, group.by = "cellType1", label = T)+NoLegend()+NoAxes()
E1618m <- NormalizeData(E1618m) %>% ScaleData() %>%FindVariableFeatures() %>%
  RunPCA(verbose = FALSE)
E1618m <- RunUMAP(E1618m, dims = 1:50, reduction.name = "umap",n.neighbors = 35L)
E1618m <- FindNeighbors(E1618m,  dims = 1:50)
E1618m <- FindClusters(object = E1618m, verbose = FALSE)
p1 <- DimPlot(E1618m, label = T)+NoAxes()+NoLegend()             
p2 <- DimPlot(E1618m, group.by = "orig.ident")+NoAxes()
p3 <- DimPlot(E1618m, group.by = "stage")+NoAxes()
p4 <- DimPlot(E1618m, group.by = "dataset")+NoAxes()
(p1+p2)/(p3+p4)
DimPlot(E1618m, label = T,)+NoAxes()+NoLegend()

#subset PC, IN and CN 
CNPCm <- subset(E1618m, idents=c("5","25","18"))
CNPCm <- NormalizeData(CNPCm) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)
CNPCm <- RunUMAP(CNPCm,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.5) %>%
  FindNeighbors( dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1.5, verbose = FALSE)
DimPlot(CNPCm, label = T)+NoAxes()+NoLegend()
FeaturePlot(CNPCm, c("Foxp2","Meis2","Pax2","Pax6"))

# Remove non-PC cells
PCm <- subset(CNPCm, idents=c(1,6,10,11,14))
PCm <- NormalizeData(PCm) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE) 
PCm <- RunUMAP(PCm,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.3)%>%
  FindNeighbors( dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 2, verbose = FALSE)
DimPlot(PCm, label = T)+NoAxes()+NoLegend()

# Load inhouse PC-only data and prepare it as the reference
PCin <- readRDS("/Users/naghamkhourifarah/Documents/PC manuscript/datta_analysis/data/PCin_new.rds")
Idents(PCin) <- "PC0.8"
DefaultAssay(PCin) <- "RNA"
PCin <- FindVariableFeatures(PCin)
PCin <- ScaleData(PCin)
PCin <- RunPCA(PCin)
PCin <- RunUMAP(PCin,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.3,reduction = "mnn")

# Use in-House PC as reference to annotate and map PCs from Taylor data
anchors <- FindTransferAnchors(reference = PCin, query = PCt,dims = 1:30, reference.reduction = "pca", features =intersect(rownames(PCin),rownames(PCt)),reduction = "pcaproject")
predictions <- TransferData(anchorset =anchors, refdata = PCin$PC0.8,
                            dims = 1:30)
query <- AddMetaData(PCt, metadata = predictions)
query <- MapQuery(anchorset = anchors, reference = PCin, query = query, refdata = list(celltype = "PC0.8"),  reduction.model = "umap")
query$predicted.id <- factor(query$predicted.id,levels = c(paste0("PC",1:11)))
query$PC0.8 <- query$predicted.id

# Use in-House PC as reference to annotate and map PCs from LaManno data
anchors1 <- FindTransferAnchors(reference = PCin, query = PCm,dims = 1:30, reference.reduction = "pca", features =intersect(rownames(PCin),rownames(PCm)),reduction = "pcaproject")
predictions <- TransferData(anchorset =anchors1, refdata = PCin$PC0.8,
                            dims = 1:30)
query1 <- AddMetaData(PCm, metadata = predictions)
query1 <- MapQuery(anchorset = anchors1, reference = PCin, query = query1, refdata = list(celltype = "PC0.8"),  reduction.model = "umap")
query1$predicted.id <- factor(query1$predicted.id,levels = c(paste0("PC",1:11)))
query1$PC0.8 <- query1$predicted.id

p1 <- DimPlot(PCin, label = T)+NoLegend()+NoAxes()
p2 <- DimPlot(query,group.by = "predicted.id",reduction = "ref.umap",label = T)+NoLegend()+NoAxes()+ggtitle("")
p3 <- DimPlot(query1,group.by = "predicted.id",reduction = "ref.umap",label = T)+NoLegend()+NoAxes()+ggtitle("")

pdf("UMAP_PC_all1.pdf", w=14,h=5)
(p1+p2+p3)
dev.off()

tab <- cbind(prop.table(table(PCin$PC0.8))*100, prop.table(table(query$predicted.id))*100, prop.table(table(query1$predicted.id))*100)
tab = as.data.frame(tab)
tab$cellType = rownames(tab)
colnames(tab) <- c("In-house","LaManno","Taylor","cellType")

df = tab %>% gather(dataset, value, c("In-house","LaManno","Taylor"))
df$cellType = factor(df$cellType, levels = gtools::mixedsort(unique(df$cellType)))

pdf("barplot_PCall.pdf", h =5, w =3)
ggplot(df, aes(x = dataset, y = value, fill = cellType))+
  geom_bar(stat="identity", color = "grey")+
  xlab("")+ylab("Percentage (%)")+
  scale_y_reverse()+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

save(PCin,query,query1, file = "data/PC_Objects.Rdata")
