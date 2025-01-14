library(Seurat)
PCin <- readRDS("data/PCin_new.rds")
cdS <- readRDS("data/E16_18_All.rds")
DefaultAssay(PCin) <- "RNA"
DimPlot(PCin)
Idents(PCin) <- "PC0.8"
PCall <- subset(cdS, idents = "PC")
PCin <- FindVariableFeatures(PCin)
PCin <- ScaleData(PCin)
PCin <- RunPCA(PCin)
PCin <- RunUMAP(PCin,dims = 1:30,return.model = T,umap.method = "uwot",n.neighbors = 35L,min.dist = 0.3,reduction = "mnn")

anchors <- FindTransferAnchors(reference = PCin, query = PCall,dims = 1:30, reference.reduction = "pca", features =intersect(rownames(PCin),rownames(PCall)))
predictions <- TransferData(anchorset =anchors, refdata = PCin$PC0.8, dims = 1:30)
query <- AddMetaData(PCall, metadata = predictions)
query <- MapQuery(anchorset = anchors, reference = PCin, query = query,refdata = list(celltype = "PC0.8"),  reduction.model = "umap")
p1 <- DimPlot(query,group.by = "predicted.id",reduction = "ref.umap",split.by = "data")+NoLegend()+NoAxes()
p2 <- DimPlot(query,group.by = "predicted.id",reduction = "ref.umap", label = T)+NoAxes()+NoLegend()
p3 <- DimPlot(cdS,label = T,label.size = 3)+NoAxes()+NoLegend()
p4 <- DimPlot(cdS,split.by = "data")+NoLegend()+NoAxes()

(p3+p2)/p4/p1

pdf("PC_all.pdf",w=8,h=12)
(p3+p2)/p4/p1
dev.off()
saveRDS(query,"query.rds")

