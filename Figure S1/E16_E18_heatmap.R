library(Seurat); library(tidyverse); library(pheatmap);library(ComplexHeatmap)
cdS_N <- readRDS("data/cdS_N.rds")
cdS_N$cellType4 <- Idents(cdS_N)
cdS_N$anno <- paste(cdS_N$stage,cdS_N$cellType4,sep = "_")
cdS_N$anno <- factor(cdS_N$anno , levels =  c(paste("E16",sort(unique(cdS_N$cellType4)),sep = "_"), c(paste("E18",sort(unique(cdS_N$cellType4)),sep = "_"))))
Idents(cdS_N) <- "cellType4"
levels(cdS_N) <- sort(levels(cdS_N))
marker_cosg <- cosg(cdS_N, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=10,remove_lowly_expressed = F)
all_markers<- marker_cosg$names %>%
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

Idents(cdS_N) <- "anno"

levels(cdS_N) <- sort(levels(cdS_N))
ave = AverageExpression(cdS_N)$RNA %>% as.matrix()
z.mat <- t(scale(t(ave), center = TRUE, scale = TRUE))
mat <- z.mat[all_markers,]

Idents(cdS_N) <- "cellType4"
levels(cdS_N) <- sort(levels(cdS_N))

p <- DimPlot(cdS_N, label = T)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata$seurat_clusters <- Idents(cdS_N)
pdata <-  pdata[order(pdata$seurat_clusters), ] # Order the plot data by group
ucols1 <- unique(pdata$colour) # Get a vector of unique colors
names(ucols1) <- levels(Idents(cdS_N))

Idents(cdS_N) <- "anno"
annotation_col <- data.frame(levels(cdS_N))
rownames(annotation_col) <- annotation_col$levels.cdS_N.
annotation_col <- annotation_col %>% separate(levels.cdS_N.,into = c("stage","cellType"),sep =  "_")
ann_colors = list(
  stage = c(E16 = "firebrick", E18="blue"),
  cellType = ucols1)
    
pheatmap(mat,cluster_rows = F,cluster_cols = F, annotation_col = annotation_col,annotation_colors = ann_colors)



paletteLength <- 50
myColor <- colorRampPalette(c("lightblue","lightblue","lightblue","white","#FFFBC8", "red","red", "firebrick"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2)))


colnames(mat) <- gsub("E16-","",colnames(mat))
colnames(mat) <- gsub("E18-","",colnames(mat))

pdf("E16_E18_markers.pdf",w=10,h=20)
pheatmap(mat,cluster_rows = F,cluster_cols = F, annotation_col = annotation_col,annotation_colors = ann_colors,show_rownames = F,color=myColor, breaks=myBreaks)
dev.off()


