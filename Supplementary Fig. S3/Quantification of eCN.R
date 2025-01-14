# =======================================================================================
# Purpose: Quantification of eCN neurons among 4 different genotypes
# Author: James Li
# Date: December 31, 2024
# 

library(Seurat);library(tidyverse)
setwd("/Users/jamesli/Documents/Manuscripts/PC_Manuscript/data_analysis/Multiplexed")

seu = readRDS("/Volumes/Jali/shared/PCpaperReviews/integrated_all.rds")

# correct the slice name of one sample
seu$ID2 = paste(seu$sample, seu$orig.ident, sep = "_")
seu$ID2 = gsub("X","", seu$ID2)
seu$ID2 = gsub("TR","P", seu$ID2)
seu$ID2 = factor(seu$ID2, levels = gtools::mixedsort(unique(seu$ID2)))

head(seu@meta.data)
table(seu$annotation)

# Remove non-cerebellar cell groups
id_remove = c("Mb.v","Mb.v6","Mb.x","Mb.v5","to.filter",
              "NPC.Mb","Mb.v1","Mb.v3","Mb.v4","Mb.v2")

Idents(seu) = "annotation"
seu_sub = subset(seu, idents = id_remove, invert = T)
table(Idents(seu_sub))

obj.list = SplitObject(seu_sub, split.by = "sample")
names <- unique(seu_sub$sample)

for (n in c("Fastigial","Interposed","Dentate","CN.gaba")){
  p = list()
  for (j in names){
    seu_tem = obj.list[[j]]
    if (n %in% levels(seu_tem)) {
      ids = WhichCells(seu_tem, idents = n)
      p[[j]] = DimPlot(seu_tem, cells.highlight = ids, reduction = "spatial", pt.size = 1, split.by = "ID2", ncol = 3, label = F,raster = F)+
        NoAxes()+NoLegend()+coord_fixed()
    }
  }
  p1 = cowplot::plot_grid(plotlist = p, ncol = 2)
  png(paste0(n,"1.png"), w=8000, h=8000)
  plot(p1)
  dev.off()
}

meta_filtered = seu_sub@meta.data
meta_filtered$annotation = droplevels(meta_filtered$annotation)

# calculate the proportion of cell groups across different samples
df = prop.table(table(meta_filtered$annotation, meta_filtered$sample), margin = 2) %>% 
  as.data.frame() %>% 
  set_names("cellType", "Sample", "Ratio")

# add a new column to the dataframe to show the percentage of the cluster in the sample
df$Percentage = df$Ratio * 100

# create a new column called genotype and add it to the dataframe
forGroup = meta_filtered[!duplicated(meta_filtered$sample),]
df$genotype = plyr::mapvalues(df$Sample, forGroup$sample,forGroup$genotype)
df$genotype = factor(df$genotype, levels = c("CTRL", "P1cKO","P2cKO","dKO"))

df_cn = df %>% 
  filter(cellType %in% c("Fastigial", "Interposed","Dentate"))

df_cn$cellType = factor(df_cn$cellType, 
                     levels = c("Fastigial", "Interposed","Dentate"))

# Create a boxplot with a jittered scatter plot overlay
p_all = ggplot(df_cn, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, ncol = 3) +
  geom_jitter(width = 0.25, size = 3, alpha = 0.6)+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "", x = "", y = "Proportion relative to total cell count (%)")
p_all


# Perform a two-way ANOVA
model <- aov(Ratio ~ cellType * genotype, data = df_cn)

# View the ANOVA table
summary(model)
# Df    Sum Sq   Mean Sq F value Pr(>F)  
# cellType           2 0.0004793 2.396e-04   3.927 0.0335 *
#   genotype           3 0.0003506 1.169e-04   1.915 0.1541  
# cellType:genotype  6 0.0000700 1.166e-05   0.191 0.9764  
# Residuals         24 0.0014645 6.102e-05                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:genotype`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("Dentate:P1cKO-Dentate:CTRL","Dentate:P2cKO-Dentate:CTRL","Dentate:dKO-Dentate:CTRL",
                                               "Interposed:P1cKO-Interposed:CTRL","Interposed:P2cKO-Interposed:CTRL","Interposed:dKO-Interposed:CTRL",
                                               "Fastigial:P1cKO-Fastigial:CTRL","Fastigial:P2cKO-Fastigial:CTRL","Fastigial:dKO-Fastigial:CTRL"
),]

filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
# [1] diff  lwr   upr   p adj

# === Repeat with eCN neuron only
celltype_selected = c("Fastigial", "Interposed","Dentate")

seu_cn = subset(seu, idents = celltype_selected)
seu_cn$annotation = droplevels(seu_cn$annotation)
# calculate the proportion of cell groups across different samples
df_cn1 = prop.table(table(seu_cn$annotation, seu_cn$sample), margin = 2) %>% 
  as.data.frame() %>% 
  set_names("cellType", "Sample", "Ratio") %>% 
  mutate(Percentage = Ratio*100)

# create a new column called genotype and add it to the dataframe
forGroup = meta_filtered[!duplicated(meta_filtered$sample),]
df_cn1$genotype = plyr::mapvalues(df_cn1$Sample, forGroup$sample,forGroup$genotype)
df_cn1$genotype = factor(df_cn1$genotype, levels = c("CTRL", "P1cKO","P2cKO","dKO"))

df_cn1$cellType = factor(df_cn1$cellType, 
                        levels = c("Fastigial", "Interposed","Dentate"))

# Create a boxplot with a jittered scatter plot overlay
p_cn = ggplot(df_cn1, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, ncol = 3) +
  geom_jitter(width = 0.25, size = 3, alpha = 0.6)+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "", x = "", y = "Proportion relative to total eCN count (%)")
p_cn

pdf("eCN_abundance.pdf", h =5, w =12)
cowplot::plot_grid(p_all, p_cn, ncol = 2)
dev.off()

# Perform a two-way ANOVA
model1 <- aov(Ratio ~ cellType * genotype, data = df_cn1)

# View the ANOVA table
summary(model1)
# Df  Sum Sq Mean Sq F value   Pr(>F)    
# cellType           2 0.07613 0.03806  12.624 0.000179 ***
#   genotype           3 0.00000 0.00000   0.000 1.000000    
# cellType:genotype  6 0.01722 0.00287   0.952 0.477460    
# Residuals         24 0.07236 0.00302                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
