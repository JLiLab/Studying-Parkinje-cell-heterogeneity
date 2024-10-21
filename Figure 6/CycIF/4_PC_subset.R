library(Seurat);library(tidyverse)
#Subset PC clusters and use scRNA data to transfer annotation
cdS <- readRDS("data/integrated_all.rds")
PC <- readRDS("data/PCi.rds")
DefaultAssay(PC) <- "RNA"
Idents(PC) <- "PC0.8"
PCsp <- subset(cdS,idents = c(paste("PC",1:11,sep = "")))
p <- DimPlot(PCsp,  pt.size = 1, label = T)+NoLegend()+NoAxes()
PCid <- CellSelector(p)
PCsp <- subset(PCsp, cells=PCid)
saveRDS(PCid, "data/PCid.rds")
allpc.id <- colnames(PCsp)
pc.filter <- setdiff(allpc.id, PCid)
cdS <-  SetIdent(cdS,cells = pc.filter, value = "to.filter" )

DimPlot(PCsp, label = T, pt.size = 1)+NoLegend()+NoAxes()
cdS$cellTypeAllMx <- Idents(cdS)
anchors <- FindTransferAnchors(reference = PC,query = PCsp, features = rownames(PCsp),reduction = "rpca")
predicted <- TransferData(anchorset = anchors,refdata = PC$PC0.8)
PCsp$predicted <- predicted$predicted.id
DimPlot(PCsp, group.by = "predicted", pt.size = 1)
Idents(PCsp) <- "predicted"
levels(PCsp) <- paste("PC",1:11,sep = "")
PCsp$predicted <- Idents(PCsp)
#PCspatial vs PCpredicted confusion Matrix
cM <- table(PCsp$predicted,PCsp$cellTypeAllM)
cM <- cM[,c(paste("PC",1:11,sep = ""))]
cM <- cM[c(paste("PC",1:11,sep = "")),]
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), clustering_method = NA,
  color = paletteContinuous("solarExtra"),treeheight_row = 0,treeheight_col = 0,
  border_color = "black",scale = "column",cluster_rows = F,cluster_cols = F)
pdf("figures/PCspatialvsPCpredicted_confusionMatrix.pdf",width = 8,height = 8)
print(p)
dev.off()
# create a column for predicted PC in the full dataset 
cdS$PC.predicted <- PCsp$predicted

# output PC-only spatial plots (I reordered the cluster to make PC2 in more visible color to show the P1-cKO phenotype)
PCsp$cellTypeAllM <- factor(PCsp$cellTypeAllM,c("PC1","PC10","PC11","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9"))
obj.list1 = SplitObject(PCsp, split.by = "sample")
names <- unique(PCsp$sample)
dir.create("figures/PCsp_cellTypeAllM")
for (j in names){
  seu_tem = obj.list1[[j]]
  p = DimPlot(seu_tem, group.by = "cellTypeAllM",reduction = "spatial", pt.size = 0.5, split.by = "ID1", ncol = 2, label = F,raster = F)+NoAxes()+NoLegend()
  
  pdf(paste0("figures/PCsp_cellTypeAllM/PCOnly_",j,".pdf"),width = 20,height =20)
  plot(p)
  dev.off()
}

# Output PC spatial plots with with other cell types in grey 
cdS$PC_cellTypeAllM <- PCsp$cellTypeAllM
Idents(cdS) <- "PC_cellTypeAllM"
p <- DimPlot(cdS, label = T)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata$seurat_clusters <- Idents(cdS)

pdata <-  pdata[order(pdata$seurat_clusters), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
names(ucols) <- levels(Idents(cdS))
ucols[12] <- "grey90"

dir.create("figures/cdS_cellTypeAllM_PC")
obj.list = SplitObject(cdS, split.by = "sample")
for (j in names){
  seu_tem = obj.list[[j]]
  p = DimPlot(seu_tem, group.by = "PC_cellTypeAllM",reduction = "spatial", pt.size = 0.5, split.by = "ID1", ncol = 2, label = F,raster = F,cols = ucols)+NoAxes()+NoLegend()
  
  pdf(paste0("figures/cdS_cellTypeAllM_PC/PCOnly_",j,".pdf"),width = 20,height =20)
  plot(p)
  dev.off()
}


Idents(cdS) <- "cellTypeAllM"
gcps <- grep("GCP", levels(cdS),value = T)
GCPC <- subset(cdS, idents=c(gcps,paste("PC",1:11,sep = "") ))
cdS$GCPC <- Idents(GCPC)
Idents(cdS) <- "GCPC"
dir.create("figures/GCPC_Only")
obj.list = SplitObject(GCPC, split.by = "sample")
for (j in names){
  seu_tem = obj.list[[j]]
  p = DimPlot(seu_tem, reduction = "spatial", pt.size = 0.5, split.by = "ID1", ncol = 2, label = F,raster = F)+NoAxes()+NoLegend()
  
  pdf(paste0("figures/GCPC_Only/GCPC_",j,".pdf"),width = 20,height =20)
  plot(p)
  dev.off()
}

Idents(cdS) <- "cellTypeAllM"
cdS$cellGroup = cdS$cellTypeAllM
meta <- cdS@meta.data
meta$cellGroup <- as.character(meta$cellGroup)
meta$cellGroup[meta$cellGroup %in% c("Dentate","Fastigial","Interposed")] <- "CN"
meta$cellGroup[grep("^GCP", meta$cellGroup)]  <- "GCP"
meta$cellGroup[meta$cellGroup %in% c(paste("GC",1:4,sep = ""))] <- "GC"
meta$cellGroup[grep("Mb", meta$cellGroup)] <- "Mb"
meta$cellGroup[grep("^PC", meta$cellGroup)] <- "PC"
table(meta$cellGroup)

cdS@meta.data <- meta
Idents(cdS) <- "cellGroup"
DimPlot(cdS,pt.size = 1,label = T)+NoLegend()
levels(cdS) <- c( "GABA.Pro","GC","Mb", "MidO","GCP","NPC" ,"to.filter","UBC", "PC","CN.gaba","GABA.Pre","IN","CN","OPC"  )
cdS$cellGroup <- Idents(cdS)
dir.create("figures/cellGroup")
obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
for (j in names){
  seu_tem = obj.list[[j]]
  p = DimPlot(seu_tem, reduction = "spatial", pt.size = 0.5, split.by = "ID1", ncol = 2, label = F,raster = F)+NoAxes()+NoLegend()
  
  pdf(paste0("figures/cellGroup/Group_",j,".pdf"),width = 20,height =20)
  plot(p)
  dev.off()
}


dir.create("figures/inspection_cellGroup")
for (n in levels(cdS)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID1", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_cellGroup/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}

head(meta)
saveRDS(cdS, "data/Integrated_all.rds")
saveRDS(PCsp, "data/PCsp.rds")
write.csv2(meta,"tables/meta_data.csv")
