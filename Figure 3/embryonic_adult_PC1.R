# integrate embryonic and adult PC scRNA-seq data

library(Seurat);library(tidyverse)

adPC <- readRDS("/Users/jamesli/Documents/Nagham/SynologyDrive/data_analysis/data/adPC.RDS")
PCin <- readRDS("/Users/jamesli/Documents/Nagham/SynologyDrive/data_analysis/data/PCin_new.rds")

DefaultAssay(adPC) <- "RNA"
ref <- adPC

ref$subcluster = gsub("Purkinje_","", ref$subcluster)
Idents(ref) <- ref$subcluster

ref <- RenameIdents(ref, "Aldoc_1"="Aldoc","Aldoc_2"="Aldoc","Aldoc_3"="Aldoc","Aldoc_4"="Aldoc","Aldoc_5"="Aldoc","Aldoc_6"="Aldoc","Aldoc_7"="Aldoc", "Anti_Aldoc_1"="Anti_Aldoc","Anti_Aldoc_2"="Anti_Aldoc")
DimPlot(ref, label = T)+NoLegend()+NoAxes()
prop.table(table(Idents(ref)))
# Aldoc Anti_Aldoc 
# 0.6268486  0.3731514 

ref$Aldoc <- Idents(ref) 
ref$Aldoc1 <- ref$subcluster 

anchors <- FindTransferAnchors(reference = ref, query = PCin, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = ref$Aldoc, weight.reduction = PCin[['mnn']],dims = 1:30)

hist(predicted.labels$prediction.score.max)
abline(v = 0.5, col = "red")

PCin <- AddMetaData(object = PCin, metadata = predicted.labels)

# Replace each label with its most likely prediction
for(i in levels(PCin)) {
  cells_to_reid <- WhichCells(PCin, idents = i)
  newid <- names(sort(table(PCin$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(PCin, cells = cells_to_reid) <- newid
}

PCin$Aldoc <- PCin$predicted.id
prop.table(table(PCin$Aldoc))
# Aldoc Anti_Aldoc 
# 0.7772512  0.2227488 

p1 <- DimPlot(PCin, group.by = "PC0.8", label = T)+NoLegend()+NoAxes()+ggtitle("PC subgroup")
p2 <- DimPlot(PCin, group.by = "Aldoc", label = T)+NoLegend()+NoAxes()+ggtitle("adult_Aldoc")

df = prop.table(table(PCin$PC0.8, PCin$Aldoc), margin = 1) %>% as.data.frame() %>% mutate(Freq = Freq*100)
colnames(df) = c("Embryonic","Adult","Percentage")

p3= ggplot(df, aes(x = Embryonic, y = Percentage, group = Adult, fill = Adult))+
  geom_bar(stat = "identity")+
  xlab("")+
  cowplot::theme_cowplot()+
  theme(legend.position = "top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("figures/embryonic_adult_PC.pdf", w =13, h=5)
cowplot::plot_grid(p1,p2,p3, rel_widths = c(1,1,0.6), ncol = 3)
dev.off()

anchors <- FindTransferAnchors(reference = ref, query = PCin, reduction = 'cca', normalization.method = "LogNormalize", dims = 1:30)
predicted.labels <- TransferData(anchorset = anchors, refdata = ref$region, weight.reduction = PCin[['mnn']],dims = 1:30)
hist(predicted.labels$prediction.score.max)
abline(v = 0.5, col = "red")

PCin <- AddMetaData(object = PCin, metadata = predicted.labels)

# Replace each label with its most likely prediction
for(i in levels(PCin)) {
  cells_to_reid <- WhichCells(PCin, idents = i)
  newid <- names(sort(table(PCin$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(PCin, cells = cells_to_reid) <- newid
}

PCin$region <- PCin$predicted.id
DimPlot(PCin, group.by = "region")+NoAxes()+ggtitle("adult_region")

df1 = prop.table(table(PCin$PC0.8, PCin$region), margin = 1) %>% as.data.frame() %>% mutate(Freq = Freq*100)
colnames(df1) = c("Embryonic","Region","Percentage")

p4= ggplot(df1, aes(x = Embryonic, y = Percentage, group = Region, fill = Region))+
  geom_bar(stat = "identity")+
  xlab("")+
  cowplot::theme_cowplot()+
  theme(legend.position = "top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

predicted.labels <- TransferData(anchorset = anchors, refdata = ref$Aldoc1, weight.reduction = PCin[['mnn']],dims = 1:30)
hist(predicted.labels$prediction.score.max)
abline(v = 0.5, col = "red")

PCin <- AddMetaData(object = PCin, metadata = predicted.labels)

# Replace each label with its most likely prediction
for(i in levels(PCin)) {
  cells_to_reid <- WhichCells(PCin, idents = i)
  newid <- names(sort(table(PCin$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(PCin, cells = cells_to_reid) <- newid
}

PCin$Aldoc1 <- PCin$predicted.id
DimPlot(PCin, group.by = "Aldoc1")+NoAxes()+ggtitle("subclusters")

df2 = prop.table(table(PCin$PC0.8, PCin$Aldoc1), margin = 1) %>% as.data.frame() %>% mutate(Freq = Freq*100)
colnames(df2) = c("Embryonic","Subcluster","Percentage")

p5= ggplot(df2, aes(x = Embryonic, y = Percentage, group = Subcluster, fill = Subcluster))+
  geom_bar(stat = "identity")+
  xlab("")+
  cowplot::theme_cowplot()+
  theme(legend.position = "top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("figures/integration.pdf", h = 5, w =8)
p4+p5
dev.off()
