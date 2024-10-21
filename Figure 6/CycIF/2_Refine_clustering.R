library(Seurat);library(tidyverse)

# Refine PC annotation using subclustering
cdS <- readRDS("../cdS_No31.rds")
# correct the slice name of one sample
cdS$ID2 = paste(cdS$orig.ident, cdS$sample, sep = "_")
cdS$ID2 = gsub("X","", cdS$ID2)
cdS$ID2 = gsub("TR","P", cdS$ID2)

Idents(cdS) <- "res1"
cdS <- RenameIdents(cdS, c("0"="GC.NPC","1"="IN","2"="PC2","3"="GC2","4"="CN.gaba","5"="PC3","6"="PC7","7"="PC1","24"="OPC","8"="Mb.v","9"="PC8",
                           "10"="GCP4","11"="Dentate.Mb","12"="Fastigial","18"="Interposed","13"="GCP1","14"="PC6","15"="UBC", "16"="NPC.Mb",
                           "17"="PC9","19"="GCP2.3", "20"="GABA.Pro","21"="PC4","22"="GC3","23"="PCx","25"="Isl1.Mb","26"="Mb.v3", "27"="GC4",
                           "28"="Mb.x","29"="PC4","30"="to.filter","31"="PC11","32"="PC7","33"="Mb.v1","34"="Mb.v2","35"="MidO","36"="Mb.v4",
                           "37"="to.filter","38"="to.filter","39"="to.filter","40"="to.filter","41"="to.filter","42"="to.filter"))
DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()

cdS <- FindSubCluster(cdS, cluster = "PC2", resolution = 1,graph.name = "RNA_snn",subcluster.name = "subPC2")
Idents(cdS) <- "subPC2"
obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_SUBPC2")
for (n in grep("PC2",levels(cdS),value = T)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_SUBPC2/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}
cdS <- RenameIdents(cdS,c("PC2_0"="PC2","PC2_1"="PC2","PC2_2"="PC2","PC2_3"="PC2","PC2_4"="PC2","PC2_5"="PC2","PC2_6"="PC2","PC2_7"="PC2",
                          "PC2_8"="PC2","PC2_9"="PC2","PC2_10"="PC9","PC2_11"="PC2","PC2_12"="PC2","PC2_13"="PC2","PC2_14"="PC11","PC2_15"="PC2",
                          "PC2_16"="PC11","PC2_17"="PC11","PC2_18"="PC11","PC2_19"="PC1"))

cdS <- FindSubCluster(cdS, cluster = "PCx", resolution = 0.1,graph.name = "RNA_snn",subcluster.name = "subPCx")
Idents(cdS) <- "subPCx"
cdS <- RenameIdents(cdS, c("PCx_0"="GABA.Pre","PCx_1"="Mb.v5","PCx_2"="Mb.v5","PCx_3"="Mb.v5","PCx_4"="GABA.Pre"))
cdS <- FindSubCluster(cdS, cluster = "PC8", resolution = 1,graph.name = "RNA_snn",subcluster.name = "subPC8")
Idents(cdS) <- "subPC8"

obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_SUBPC8")
for (n in grep("PC8",levels(cdS),value = T)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_SUBPC8/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}

cdS <- RenameIdents(cdS, c("PC8_0"="PC11","PC8_1"="PC8","PC8_2"="PC8","PC8_3"="PC5","PC8_4"="PC8","PC8_5"="PC8","PC8_6"="PC8",
                           "PC8_7"="PC5","PC8_8"="PC5","PC8_9"="PC8","PC8_10"="PC8","PC8_11"="PC5","PC8_12"="PC8","PC8_13"="PC11",
                           "PC8_14"="PC11", "PC8_15"="PC5","PC8_16"="PC8","PC8_17"="PC9","PC8_18"="PC8"))

cdS <- FindSubCluster(cdS, cluster = "PC3", resolution =1,graph.name = "RNA_snn",subcluster.name = "subPC3")
Idents(cdS) <- "subPC3"
VlnPlot(cdS, c("Calb2","Pcdh10","Nr2f2","Plcb4","Nrgn"),stack = T,pt.size = 0,flip = TRUE,sort = T,idents = c(grep("PC3",levels(cdS),value = T)))

obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_SUBPC3")
for (n in grep("PC3",levels(cdS),value = T)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_SUBPC3/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}

cdS <- RenameIdents(cdS, c("PC3_0"="PC10","PC3_1"="PC3","PC3_2"="PC3","PC3_3"="PC3","PC3_4"="PC3","PC3_5"="PC3","PC3_6"="PC9","PC3_7"="PC3","PC3_8"="PC9","PC3_9"="PC3",
                           "PC3_10"="PC9","PC3_11"="PC3","PC3_12"="PC3","PC3_13"="PC9","PC3_14"="PC3","PC3_15"="PC10","PC3_16"="PC9","PC3_17"="PC3","PC3_18"="PC3"))

cdS <- FindSubCluster(cdS, cluster = "PC9", resolution = 0.5,graph.name = "RNA_snn",subcluster.name = "subPC9")
Idents(cdS) <- "subPC9"

obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_SUBPC9")
for (n in grep("PC9",levels(cdS),value = T)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_SUBPC9/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}

cdS <- RenameIdents(cdS, c("PC9_0"="PC9","PC9_1"="PC9","PC9_2"="PC7","PC9_3"="CN.gaba","PC9_4"="PC9","PC9_5"="PC2",
                           "PC9_6"="PC9","PC9_7"="PC9","PC9_8"="PC3","PC9_9"="PC3" ))


cdS$cellTypeAllM2 <- Idents(cdS)
DimPlot(cdS,  pt.size = 1, label = T)+NoLegend()+NoAxes()

#Refine GCP clusters
cdS <- FindSubCluster(cdS, cluster = "GCP2.3", resolution = 0.5,graph.name = "RNA_snn",subcluster.name = "subGCP2.3")
Idents(cdS) <- "subGCP2.3"
cdS <- RenameIdents(cdS, c("GCP2.3_0"="GCP.Tlx3","GCP2.3_1"="GCP.Tlx3","GCP2.3_2"="GCP.Tlx3","GCP2.3_3"="GCP.Tlx3","GCP2.3_4"="GCP.Isl1","GCP2.3_5"="GCP.Tlx3","GCP2.3_6"="GCP.Isl1",
                           "GCP2.3_7"="GCP.Tlx3","GCP2.3_8"="GCP.Tlx3","GCP2.3_9"="GCP.Tlx3","GCP2.3_10"="GCP.Tlx3","GCP2.3_11"="GCP.Tlx3"))

cdS <- FindSubCluster(cdS, cluster = "GCP4", resolution = 0.2,graph.name = "RNA_snn",subcluster.name = "subGCP4")
Idents(cdS) <- "subGCP4"
cdS <- RenameIdents(cdS, c("GCP.Tlx3" ="GCP2.3","GCP4_0"="GCP4","GCP4_1"="UBC","GCP4_2"="GCP4","GCP4_3"="GCP4","GCP4_4"="GCP4"))
DimPlot(cdS, label = T, pt.size = 1)+NoLegend()+NoAxes()
cdS$cellTypeAllM1 <- Idents(cdS)

#refine clustering of other cell types
cdS <- FindSubCluster(cdS, cluster = "Dentate.Mb", resolution = 0.5,graph.name = "RNA_snn",subcluster.name = "Dent")
Idents(cdS) <- "Dent"
cdS <- RenameIdents(cdS, c("Dentate.Mb_0" ="Dentate","Dentate.Mb_1" ="Dentate","Dentate.Mb_2" ="Dentate","Dentate.Mb_3" ="Dentate",
                           "Dentate.Mb_7" ="Dentate","Dentate.Mb_4" ="Mb.v6","Dentate.Mb_5" ="Mb.v6","Dentate.Mb_6" ="Mb.v6",
                           "Dentate.Mb_8" ="Mb.v6","Dentate.Mb_9" ="Mb.v6","Dentate.Mb_10" ="Mb.v6","Dentate.Mb_11" ="Mb.v6","Dentate.Mb_12" ="Mb.v6","Dentate.Mb_13" ="Mb.v6"))

cdS <- FindSubCluster(cdS, cluster = "GC.NPC", resolution = 0.5,graph.name = "RNA_snn",subcluster.name = "GC")
Idents(cdS) <- "GC"
cdS <- RenameIdents(cdS, c("GC.NPC_0" ="GC1","GC.NPC_1" ="GC1","GC.NPC_2" ="NPC","GC.NPC_9" ="NPC","GC.NPC_3" ="GC1",
                           "GC.NPC_4" ="GC1","GC.NPC_5" ="GC1","GC.NPC_6" ="GC1","GC.NPC_7" ="GC1","GC.NPC_8" ="GC1",
                           "GC.NPC_10" ="GC1","GC.NPC_1" ="GC1","GC.NPC_11" ="GC1","GC.NPC_12" ="GC1"))

cdS$cellTypeAllM1 <- Idents(cdS)

Idents(cdS) <- "celltypeM"
id <- WhichCells(cdS, idents = "PC11")
Idents(cdS) <- "cellTypeAllM1"
cdS <- SetIdent(cdS, cells = id, value = 'PC11')
cdS$cellTypeAllM <- Idents(cdS)
#inspect clusters 
Idents(cdS) <- "cellTypeAllM"
obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_cellTypeAllM")
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
  png(paste0("figures/inspection_cellTypeAllM/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}



cdS <- FindSubCluster(cdS, cluster = "PC1", resolution = 0.5,graph.name = "RNA_snn",subcluster.name = "subPC1")
Idents(cdS) <- "subPC1"
levels(cdS)
obj.list = SplitObject(cdS, split.by = "sample")
names <- unique(cdS$sample)
dir.create("figures/inspection_SUBPC1")
for (n in grep("PC1",levels(cdS),value = T)){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+NoAxes()+NoLegend()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 1)
  png(paste0("figures/inspection_SUBPC1/cluster_",n,".png"), w=3000, h=10000)
  plot(p1)
  dev.off()
}
cdS <- RenameIdents(cdS, c("PC1_0"="PC1","PC1_1"="PC1","PC1_2"="PC1","PC1_3"="PC1","PC1_4"="PC1","PC1_5"="PC1",
                           "PC1_6"="PC1","PC1_7"="PC1","PC1_8"="PC11","PC1_9"="PC1","PC1_10"="PC1","PC1_11"="PC1","PC1_12"="PC1"))


cdS$cellTypeAllM <- Idents(cdS)

# Clean up meta data
meta <- cdS@meta.data
head(meta)
dim(meta)
meta = meta[,c(1:18,61,63,64,68,72,73,74,83,86)]
head(meta)
cdS@meta.data <- meta
saveRDS(cdS, "data/integrated_all.rds")
saveRDS(meta, "data/meta__filtered_new.rds")
