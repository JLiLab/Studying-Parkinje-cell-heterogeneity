library(Seurat)
library(tidyverse)

# require(biomaRt)
# chicken <- useMart('ensembl', dataset = 'ggallus_gene_ensembl', host="https://useast.ensembl.org")
# mouse <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl', host="https://useast.ensembl.org")
# annot_table <- getLDS(
#   mart = mouse,
#   attributes = c('ensembl_gene_id','mgi_symbol','external_gene_name','chromosome_name'),
#   martL = chicken,
#   attributesL = c('ensembl_gene_id','external_gene_name','chromosome_name','gene_biotype'))


pc_mouse <- readRDS("/Users/jamesli/Library/CloudStorage/OneDrive-UConnHealthCenter/Manuscript/Purkinje\ cells/addition/Chick/PCi.rds")
Idents(pc_mouse) = pc_mouse$PC0.8
DimPlot(pc_mouse, label=T)+NoLegend()+NoAxes()

seu_chick <-readRDS( "data/seurat_filtered.rds")
seu_human = readRDS("/Volumes/SSD_JamesLi/Experiments1/scRNAseq/human/Aldinger/all/data/mouse_orth.rds")

Idents(seu_chick) = "seurat_clusters"
pc_chick <- subset(seu_chick, idents = c("1","3","4","9","12","13"))

# quick conversion of chick gene symbol to mouse
featuresM <- rownames(pc_mouse)
featuresMC <- toupper(featuresM)
tab <- as.data.frame(featuresM)
head(tab)
tab$featuresMC <- featuresMC

common <- intersect(rownames(pc_chick),featuresMC)
length(common)
pc_chick <- pc_chick[common,] 

featuresMs <- tab$featuresM[match(rownames(pc_chick), tab$featuresMC)]
rownames(pc_chick) <- featuresMs

anchors <- FindTransferAnchors(reference = pc_mouse, query = pc_chick, dims = 1:30, reference.reduction = "pca", features =NULL,reduction = "pcaproject")
predictions <- TransferData(anchorset =anchors, refdata = pc_mouse$PC0.8, dims = 1:30)
query <- AddMetaData(pc_chick, metadata = predictions)
query <- MapQuery(anchorset = anchors, reference = pc_mouse, query = query,refdata = list(celltype = "PC0.8"),  reduction.model = "umap")

all(colnames(pc_chick) == colnames(query))
pc_chick$predicted.celltype = query$predicted.celltype

meta1 = data.frame(species = "chick", cellType = query$predicted.id, dataset = query$sample)
meta2 = data.frame(species = "mouse", cellType = pc_mouse$PC0.8, dataset = pc_mouse$dataset)
meta = rbind(meta1, meta2)

meta_pca = meta %>% group_by(species,dataset,cellType) %>% summarise(count = n())
meta_pca= meta_pca %>% spread(key = cellType, value = count)
# head(meta_pca)

meta_pca$PC11 <- as.numeric(meta_pca$PC11)
meta_pca$PC8 <- as.numeric(meta_pca$PC8)
meta_pca[is.na(meta_pca)] <- 0.5


meta$cellType = factor(meta$cellType, levels = paste0("PC", 1:11))

# calculate the proportion of cell clusters across different samples
df = prop.table(table(meta$cellType, meta$dataset), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(df) <- c("cellType", "dataset", "Ratio")

# add a new column to the dataframe to show the percentage of the cluster in the sample
df$Percentage = df$Ratio * 100

forGroup = meta[!duplicated(meta$dataset),]
# create a new column called genotype and add it to the dataframe
df$species = plyr::mapvalues(df$dataset, forGroup$dataset,forGroup$species)
head(df)

p = ggplot(df, aes(x = species, y = Percentage, fill = species)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = dataset, size = 1), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free") +
  # Set the color palette and legend
  scale_fill_discrete(name = "cellType") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # Add axis and title labels
  labs(title = "PC subtype Composition", x = "Species", y = "Percentage")
p

seu_human <- readRDS("/Volumes/SSD_JamesLi/Experiments1/scRNAseq/human/Aldinger/all/data/mouse_orth.rds")

# I ran into an issue to subset data in Seurat version 5. Therefore, I re-created a new Seurat object to circumvent the problem.
seu_human = CreateSeuratObject(counts = LayerData(seu_human, layer = "counts"), meta.data = seu_human@meta.data)

cell_ids = seu_human@meta.data %>% filter(fig_cell_type == "H-PC") %>% rownames(); length(cell_ids)
pc_human <- subset(seu_human, cells = cell_ids)
pc_human

features = intersect(rownames(pc_mouse), rownames(pc_human))
length(features)

# # I ran into errors with "LogNormalize" normalization and I switched to SCTransformation
# anchorsMH <- FindTransferAnchors(reference = pc_mouse, query = seu_human_subset, dims = 1:30, reference.reduction = "pca", features =NULL,reduction = "pcaproject")
# predictionsMH <- TransferData(anchorset = anchorsMH, refdata = pc_mouse$PC0.8, dims = 1:30)
# queryH <- AddMetaData(pc_human, metadata = predictionsH)

DefaultAssay(pc_mouse) <- "RNA"
pc_mouse <- SCTransform(pc_mouse, verbose = F) %>% 
  RunPCA(verbose = F)

pc_human <- SCTransform(pc_human, verbose = F) %>% 
  RunPCA(verbose = F) 

anchors <- FindTransferAnchors(
  reference = pc_mouse, 
  reference.assay = "SCT",
  query = pc_human,  
  dims = 1:30, 
  reference.reduction = "pca", 
  query.assay = "SCT",
  normalization.method = "SCT",
  recompute.residuals = FALSE
)

pc_human_predicted <- MapQuery(
  anchorset = anchors,
  query = pc_human,
  reference = pc_mouse,
  refdata = "PC0.8",
  reference.reduction = "pca"
)

all(colnames(pc_human) == colnames(pc_human_predicted))
pc_human$predicted.celltype = pc_human_predicted$predicted.id

meta3 = data.frame(species = "human", cellType = pc_human_predicted$predicted.id, dataset = pc_human_predicted$age)
meta_all = rbind(meta, meta3)

meta_all$cellType = factor(meta_all$cellType, levels = paste0("PC", 1:11))

meta_all1 = meta_all %>% filter(!dataset %in% c("14 PCW","17 PCW","16 PCW","18 PCW","20 PCW", "Cb_E6"))
meta_all1$cellType = factor(meta_all1$cellType, levels = paste0("PC", 1:11))

# save objects
save(pc_human, pc_chick, meta_all, file = "data/PC_integration.rdata")

# calculate the proportion of cell clusters across different samples
df = prop.table(table(meta_all1$cellType, meta_all1$dataset), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(df) <- c("cellType", "stage", "Ratio")

# add a new column to the dataframe to show the percentage of the cluster in the sample
df$Percentage = df$Ratio * 100

forGroup = meta_all1[!duplicated(meta_all1$dataset),]
# create a new column called genotype and add it to the dataframe
df$species = plyr::mapvalues(df$stage, forGroup$dataset,forGroup$species)
df$species = factor(df$species, levels = c("chick", "mouse", "human"))
head(df)

p = ggplot(df, aes(x = species, y = Percentage, fill = species)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = stage, size = 1), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # Add axis and title labels
  labs(title = "PC subtype Composition", x = "Species", y = "Percentage of total PCs")

pdf("figures/PC_species_test.pdf", w = 9, h =6)
p
dev.off()

df_pc1 = df %>% filter(cellType == "PC1")
p_pc1 = ggplot(df_pc1, aes(x = species, y = Percentage, fill = species)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = stage, size = 1), position = position_jitter(width = 0.3, height = 0)) +
  cowplot::theme_cowplot()+
  xlab("")+ylab("PC1 abundance (%)")+
  theme(legend.position = "none")

pdf("figures/PC1_abundance.pdf", w = 3, h = 4)
p_pc1
dev.off()


# Perform a two-way ANOVA
df$species = factor(df$species, levels = c("chick", "human","mouse"))
model <- aov(Ratio ~ cellType * species, data = df)

# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:species`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[grep("PC[0-11]+:mouse - PC[0-11]+:human", rownames(specific_comparisons)), ]

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("PC1:mouse-PC1:human","PC2:mouse-PC2:human",
                                               "PC3:mouse-PC3:human","PC4:mouse-PC4:human",
                                               "PC5:mouse-PC5:human","PC6:mouse-PC6:human",
                                               "PC7:mouse-PC7:human","PC8:mouse-PC8:human",
                                               "PC9:mouse-PC9:human","PC10:mouse-PC10:human",
                                               "PC11:mouse-PC11:human", 
                                               "PC1:mouse-PC1:chick","PC2:mouse-PC2:chick",
                                               "PC3:mouse-PC3:chick","PC4:mouse-PC4:chick",
                                               "PC5:mouse-PC5:chick","PC6:mouse-PC6:chick",
                                               "PC7:mouse-PC7:chick","PC8:mouse-PC8:chick",
                                               "PC9:mouse-PC9:chick","PC10:mouse-PC10:chick",
                                               "PC11:mouse-PC11:chick",
                                               "PC1:human-PC1:chick","PC2:human-PC2:chick",
                                               "PC3:human-PC3:chick","PC4:human-PC4:chick",
                                               "PC5:human-PC5:chick","PC6:human-PC6:chick",
                                               "PC7:human-PC7:chick","PC8:human-PC8:chick",
                                               "PC9:human-PC9:chick","PC10:human-PC10:chick",
                                               "PC11:human-PC11:chick"), ]
filtered_comparisons=as.data.frame(filtered_comparisons)

# View the filtered comparisons
print(filtered_comparisons[filtered_comparisons$'p adj'<0.05,])
#                           diff          lwr         upr        p adj
# PC1:mouse-PC1:human -0.2015664 -0.318547007 -0.08458583 9.982751e-07
# PC1:mouse-PC1:chick  0.1273492  0.002291703  0.25240672 4.041512e-02
# PC2:mouse-PC2:chick  0.1441457  0.019088215  0.26920323 7.361785e-03
# PC4:mouse-PC4:chick -0.4444307 -0.569488181 -0.31937317 0.000000e+00
# PC1:human-PC1:chick  0.3289156  0.211935044  0.44589622 0.000000e+00
# PC2:human-PC2:chick  0.1282250  0.011244392  0.24520557 1.550552e-02
# PC4:human-PC4:chick -0.3515993 -0.468579876 -0.23461870 0.000000e+00