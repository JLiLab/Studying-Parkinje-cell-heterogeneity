library(Seurat);library(tidyverse)

# reload data
cdS <- readRDS("data/integrated_all.rds")
cdS = readRDS("/Volumes/SSD_JamesLi/Experiments1/Multiplexed/Integrated_all/data/integrated_all.rds")

metaDD <- cdS@meta.data
# remove the posterior-most sections and cell clusters of imaging artifacts 
# we also removed midbrain cells because their number changed 
# (similar results obtained without filtering)
metaDD$cellTypeAllM = as.character(metaDD$cellTypeAllM)
meta_filtered = metaDD %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(!cellTypeAllM %in% c("to.filter","MidO","NPC.Mb")) %>% 
  filter(!startsWith(cellTypeAllM, "Mb."))

# groups cell types into major cell groups
table(meta_filtered$cellTypeAllM)

meta_filtered$cellGroup = meta_filtered$cellTypeAllM
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("Dentate","Fastigial","Interposed")] <- "CN"
meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GCP"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("GC1","GC2","GC3","GC4")] <- "GC"
#meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("IN","GABA.Pre")] <- "IN"
# meta_filtered$cellGroup[grep("Mb", meta_filtered$cellGroup)] <- "Mb"
meta_filtered$cellGroup[grep("^PC", meta_filtered$cellGroup)] <- "PC"
table(meta_filtered$cellGroup)

meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GCP/GC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("GC1","GC2","GC3","GC4")] <- "GCP/GC"
meta_filtered$cellGroup[meta_filtered$cellGroup == "UBC"]  <- "GCP/GC"

# calculate the proportion of cell groups across different samples
df = prop.table(table(meta_filtered$cellGroup, meta_filtered$sample), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(df) <- c("cellType", "Sample", "Ratio")

# add a new column to the dataframe to show the percentage of the cluster in the sample
df$Percentage = df$Ratio * 100

# create a new column called genotype and add it to the dataframe
forGroup = meta_filtered[!duplicated(meta_filtered$sample),]
df$genotype = plyr::mapvalues(df$Sample, forGroup$sample,forGroup$genotype)
df$genotype = factor(df$genotype, levels = c("CTRL", "P1cKO","P2cKO","dKO"))
df$cellType = factor(df$cellType, 
                     levels = c("PC","IN","CN.gaba","GABA.Pre",
                                "GC","GCP","CN","UBC",
                                "GABA.Pro","NPC","OPC","Isl1.Mb"))

df$cellType = factor(df$cellType, 
                     levels = c("PC","IN","CN.gaba","GABA.Pre",
                                "GC","GCP","CN","UBC",
                                "GABA.Pro","NPC","OPC","Isl1.Mb"))

# Create a boxplot with a jittered scatter plot overlay
p_all = ggplot(df, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, scales = "free", ncol = 4) +
  geom_jitter(width = 0.25, size = 3, alpha = 0.6)+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "Relative abundance of major cellTypes", x = "cellType", y = "Percentage")
p_all

pdf("CellTypes_abundance_all.pdf", h =8, w =10)
p_all
dev.off()

# Perform a two-way ANOVA
model <- aov(Ratio ~ cellType * genotype, data = df)

# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:genotype`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("GABA.Pro:P1cKO-GABA.Pro:CTRL","GABA.Pro:P2cKO-GABA.Pro:CTRL","GABA.Pro:dKO-GABA.Pro:CTRL",
                                               "PC:P1cKO-PC:CTRL","PC:P2cKO-PC:CTRL","PC:dKO-PC:CTRL",
                                               "IN:P1cKO-IN:CTRL","IN:P2cKO-IN:CTRL","IN:dKO-IN:CTRL",
                                               "CN.gaba:P1cKO-CN.gaba:CTRL","CN.gaba:P2cKO-CN.gaba:CTRL","CN.gaba:dKO-CN.gaba:CTRL",
                                               "CN:P1cKO-CN:CTRL","CN:P2cKO-CN:CTRL","CN:dKO-CN:CTRL",
                                               "GCP:P1cKO-GCP:CTRL","GCP:P2cKO-GCP:CTRL","GCP:dKO-GCP:CTRL",
                                               "GABA.Pre:P1cKO-GABA.Pre:CTRL","GABA.Pre:P2cKO-GABA.Pre:CTRL","GABA.Pre:dKO-GABA.Pre:CTRL",
                                               "UBC:P1cKO-UBC:CTRL","UBC:P2cKO-UBC:CTRL","UBC:dKO-UBC:CTRL",
                                               "NPC:P1cKO-NPC:CTRL","NPC:P2cKO-NPC:CTRL","NPC:dKO-NPC:CTRL",
                                               "GC:P1cKO-GC:CTRL","GC:P2cKO-GC:CTRL","GC:dKO-GC:CTRL"
),]

filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
# [1] diff  lwr   upr   p adj
# <0 rows> (or 0-length row.names)


# Examine the relative abundance of PC subtypes
meta_pc = metaDD %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(startsWith(cellTypeAllM, "PC"))

tab = prop.table(table(meta_pc$cellTypeAllM, meta_pc$sample), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(tab) <- c("cellType", "sample", "Ratio")
tab$cellType = factor(tab$cellType, levels = paste0("PC",1:11))

# add a new column to the dataframe to show the percentage of the cluster in the sample
tab$Percentage = tab$Ratio * 100

forGroup = meta_pc[!duplicated(meta_pc$sample),]
forGroup$genotype1 = forGroup$genotype
forGroup$genotype1[forGroup$genotype1 %in% c("P2cKO", "dKO")] = "P2.deficient"

# create a new column called genotype and add it to the dataframe
tab$genotype = plyr::mapvalues(tab$sample, forGroup$sample,forGroup$genotype)
tab$genotype1 = plyr::mapvalues(tab$sample, forGroup$sample,forGroup$genotype1)
tab$genotype = factor(tab$genotype, levels = c("CTRL","P1cKO","P2cKO","dKO"))
tab$genotype1 = factor(tab$genotype1, levels = c("CTRL","P1cKO","P2.deficient"))

p_pc = ggplot(tab, aes(x = genotype1, y = Percentage, fill = genotype1)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = genotype, size = 0.5), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free") +
  # Set the color palette and legend
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "Relative abundance of PC subtypes", x = "cellType", y = "Percentage")
p_pc

pdf("CellTypes_abundance_PC.pdf", h =8, w =10)
p_pc
dev.off()

# Perform a two-way ANOVA
model <- aov(Ratio ~ cellType * genotype, data = tab)

# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:genotype`


# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("PC1:P1cKO-PC1:CTRL","PC1:P2cKO-PC1:CTRL","PC1:dKO-PC1:CTRL",
                                               "PC2:P1cKO-PC2:CTRL","PC2:P2cKO-PC2:CTRL","PC2:dKO-PC2:CTRL",
                                               "PC3:P1cKO-PC3:CTRL","PC3:P2cKO-PC3:CTRL","PC3:dKO-PC3:CTRL",
                                               "PC4:P1cKO-PC4:CTRL","PC4:P2cKO-PC4:CTRL","PC4:dKO-PC4:CTRL",
                                               "PC5:P1cKO-PC5:CTRL","PC5:P2cKO-PC5:CTRL","PC5:dKO-PC5:CTRL",
                                               "PC6:P1cKO-PC6:CTRL","PC6:P2cKO-PC6:CTRL","PC6:dKO-PC6:CTRL",
                                               "PC7:P1cKO-PC7:CTRL","PC7:P2cKO-PC7:CTRL","PC7:dKO-PC7:CTRL",
                                               "PC8:P1cKO-PC8:CTRL","PC8:P2cKO-PC8:CTRL","PC8:dKO-PC8:CTRL",
                                               "PC11:P1cKO-PC11:CTRL","PC11:P2cKO-PC11:CTRL","PC11:dKO-PC11:CTRL",
                                               "PC9:P1cKO-PC9:CTRL","PC9:P2cKO-PC9:CTRL","PC9:dKO-PC9:CTRL"
),]

filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
#                           diff          lwr          upr        p adj
# PC1:P2cKO-PC1:CTRL -0.07682066 -0.134535511 -0.019105808 3.954792e-04
# PC1:dKO-PC1:CTRL   -0.10873799 -0.166452841 -0.051023138 1.934569e-08
# PC7:P2cKO-PC7:CTRL -0.07555849 -0.133273343 -0.017843640 5.652943e-04
# PC7:dKO-PC7:CTRL   -0.06124622 -0.118961076 -0.003531373 2.285155e-02
# PC11:dKO-PC11:CTRL  0.06098195  0.003267099  0.118696802 2.428300e-02

# Perform a two-way ANOVA
model <- aov(Ratio ~ cellType * genotype1, data = tab)
# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:genotype1`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("PC1:P1cKO-PC1:CTRL","PC1:P2.deficient-PC1:CTRL",
                                               "PC2:P1cKO-PC2:CTRL","PC2:P2.deficient-PC2:CTRL",
                                               "PC3:P1cKO-PC3:CTRL","PC3:P2.deficient-PC3:CTRL",
                                               "PC4:P1cKO-PC4:CTRL","PC4:P2.deficient-PC4:CTRL",
                                               "PC5:P1cKO-PC5:CTRL","PC5:P2.deficient-PC5:CTRL",
                                               "PC6:P1cKO-PC6:CTRL","PC6:P2.deficient-PC6:CTRL",
                                               "PC7:P1cKO-PC7:CTRL","PC7:P2.deficient-PC7:CTRL",
                                               "PC8:P1cKO-PC8:CTRL","PC8:P2.deficient-PC8:CTRL",
                                               "PC11:P1cKO-PC11:CTRL","PC11:P2.deficient-PC11:CTRL",
                                               "PC9:P1cKO-PC9:CTRL","PC9:P2.deficient-PC9:CTRL"
),]

filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
#                                    diff          lwr         upr        p adj
# PC1:P2.deficient-PC1:CTRL   -0.09277932 -0.145092272 -0.04046638 2.129528e-07
# PC7:P2.deficient-PC7:CTRL   -0.06840236 -0.120715305 -0.01608941 6.822051e-04
# PC11:P2.deficient-PC11:CTRL  0.05388699  0.001574043  0.10619994 3.479866e-02

