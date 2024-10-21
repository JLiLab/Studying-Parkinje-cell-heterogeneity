library(Seurat);library(tidyverse)
#Subset PC clusters and use scRNA data to transfer annotation
cdS <- readRDS("data/integrated_all.rds")

Idents(cdS) <- "cellTypeAllM"
obj.list = SplitObject(cdS, split.by = "sample")

# print CTRL sections
for (n in levels(Idents(cdS))){
  temp = obj.list[["X43430"]]
  cell_id <- WhichCells(temp, idents = n)
  plot1 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X46689"]]
  cell_id <- WhichCells(temp, idents = n)
  plot2 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X53168"]]
  cell_id <- WhichCells(temp, idents = n)
  plot3 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  plot = plot1 + plot2 + plot3
  
  png(paste0("figures/inspection_cellTypeAllM/",n,"_CTRL.png"), w=3000, h=2800)
  print(plot)
  dev.off()
}

# print P1cKO sections
for (n in levels(Idents(cdS))){
  temp = obj.list[["X43429"]]
  cell_id <- WhichCells(temp, idents = n)
  plot1 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X45820"]]
  cell_id <- WhichCells(temp, idents = n)
  plot2 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X45824"]]
  cell_id <- WhichCells(temp, idents = n)
  plot3 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  plot = plot1 + plot2 + plot3
  
  png(paste0("figures/inspection_cellTypeAllM/",n,"_P1cKO.png"), w=3000, h=2800)
  print(plot)
  dev.off()
}

cellTypes = setdiff(levels(Idents(cdS)), "Mb.v4")

# print P2cKO sections
for (n in cellTypes){
  temp = obj.list[["X45650"]]
  cell_id <- WhichCells(temp, idents = n)
  plot1 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X49530"]]
  cell_id <- WhichCells(temp, idents = n)
  plot2 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X53218"]]
  cell_id <- WhichCells(temp, idents = n)
  plot3 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  plot = plot1 + plot2 + plot3
  
  png(paste0("figures/inspection_cellTypeAllM/",n,"_P2cKO.png"), w=3000, h=2800)
  print(plot)
  dev.off()
}

cellTypes = setdiff(levels(Idents(cdS)), c("Mb.v4", "Mb.v2"))
# print P12dKO sections
for (n in cellTypes){
  temp = obj.list[["X44675"]]
  cell_id <- WhichCells(temp, idents = n)
  plot1 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X44677"]]
  cell_id <- WhichCells(temp, idents = n)
  plot2 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  temp = obj.list[["X53166"]]
  cell_id <- WhichCells(temp, idents = n)
  plot3 = DimPlot(temp, reduction = "spatial", cells.highlight = cell_id, pt.size = 0.8, split.by = "ID1", ncol = 1, label = F,raster = F)+
    NoAxes()+NoLegend()+ coord_fixed()
  
  plot = plot1 + plot2 + plot3
  
  png(paste0("figures/inspection_cellTypeAllM/",n,"_P12dKO.png"), w=3000, h=2800)
  print(plot)
  dev.off()
}


