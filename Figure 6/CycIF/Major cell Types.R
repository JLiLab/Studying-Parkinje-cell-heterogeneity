

library(Seurat);library(tidyverse)

# reload data
# cdS <- readRDS("data/integrated_all.rds")
cdS = readRDS("/Volumes/SSD_JamesLi/Experiments1/Multiplexed/Integrated_all/data/integrated_all.rds")

metaDD <- cdS@meta.data
# remove the posterior-most sections and cell clusters of imaging artifacts 
# we also removed midbrain cells because their number changed 
# (similar results obtained without filtering)
metaDD$cellTypeAllM = as.character(metaDD$cellTypeAllM)
meta_filtered = metaDD %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(!cellTypeAllM %in% c("to.filter","MidO","NPC.Mb","Isl1.Mb")) %>% 
  filter(!startsWith(cellTypeAllM, "Mb."))

# groups cell types into major cell groups
table(meta_filtered$cellTypeAllM)

meta_filtered$cellGroup = meta_filtered$cellTypeAllM
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("Dentate","Fastigial","Interposed")] <- "CN"
meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GCP"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("GC1","GC2","GC3","GC4")] <- "GC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("IN","GABA.Pre")] <- "IN"
# meta_filtered$cellGroup[grep("Mb", meta_filtered$cellGroup)] <- "Mb"
meta_filtered$cellGroup[grep("^PC", meta_filtered$cellGroup)] <- "PC"
meta_filtered$cellGroup[meta_filtered$cellGroup == "NPC"] <- "NPC/AS"

meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GCP/GC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("GC1","GC2","GC3","GC4","UBC","GC")] <- "GCP/GC"
table(meta_filtered$cellGroup)

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
                     levels = c("PC","IN","CN.gaba",
                                "GCP/GC","CN","GABA.Pro",
                                "NPC/AS","OPC"))

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

pdf("CellGroup_abundance.pdf", h =5, w =10)
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
                                               "GCP/GC:P1cKO-GCP/GC:CTRL","GCP/GC:P2cKO-GCP/GC:CTRL","GCP/GC:dKO-GCP/GC:CTRL",
                                               "NPC/AS:P1cKO-NPC/AS:CTRL","NPC/AS:P2cKO-NPC/AS:CTRL","NPC/AS:dKO-NPC/AS:CTRL",
                                               "OPC:P1cKO-OPC:CTRL","OPC:P2cKO-OPC:CTRL","OPC:dKO-OPC:CTRL"),]


filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]


# == repeat following the multiome data 
meta <- cdS@meta.data
# remove the posterior-most sections and cell clusters of imaging artifacts 
# we also removed midbrain cells because their number changed 
# (similar results obtained without filtering)
meta$cellTypeAllM = as.character(meta$cellTypeAllM)
meta_filtered = meta %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(!cellTypeAllM %in% c("to.filter","MidO","NPC.Mb","Isl1.Mb","OPC")) %>% 
  filter(!startsWith(cellTypeAllM, "Mb."))

# groups cell types into major cell groups
table(meta_filtered$cellTypeAllM)

meta_filtered$cellGroup = meta_filtered$cellTypeAllM
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("Dentate","Fastigial","Interposed")] <- "CN"
meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("GC1","GC2","GC3","GC4","UBC")] <- "GC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("IN","GABA.Pre")] <- "IN"
meta_filtered$cellGroup[grep("^PC", meta_filtered$cellGroup)] <- "PC"
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("NPC","GABA.Pro")] <- "NPC"
table(meta_filtered$cellGroup)

# calculate the proportion of cell groups across different samples
df = prop.table(table(meta_filtered$cellGroup, meta_filtered$sample), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(df) <- c("cellClass", "Sample", "Ratio")

# add a new column to the dataframe to show the percentage of the cluster in the sample
df$Percentage = df$Ratio * 100

# create a new column called genotype and add it to the dataframe
forGroup = meta_filtered[!duplicated(meta_filtered$sample),]
df$genotype = plyr::mapvalues(df$Sample, forGroup$sample,forGroup$genotype)
df$genotype = factor(df$genotype, levels = c("CTRL", "P1cKO","P2cKO","dKO"))

df$group = as.character(df$genotype)
df$group[df$group %in% c("P2cKO", "dKO")] = "P2.deficient"
df$group[df$group %in% c("CTRL", "P1cKO")] = "CTRL"
df$group = factor(df$group, levels = c("CTRL", "P2.deficient"))

df_order = df %>% 
  filter(genotype == "CTRL") %>% 
  group_by(cellClass) %>%
  summarize(mean = mean(Ratio)) %>% 
  arrange(desc(mean))
df$cellClass = factor(df$cellClass, levels = df_order$cellClass)

# Create a boxplot with a jittered scatter plot overlay
p_all = ggplot(df, aes(x = genotype, y=Percentage,fill = genotype))+
  facet_wrap(~cellClass, scales = "free", ncol = 3)+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size = 3, alpha = 0.6)+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  labs(title = "Relative abundance of major cellTypes", x = "cellType", y = "Percentage")
p_all

pdf("cellClass_abundance.pdf", h =6, w = 7)
p_all
dev.off()


# Create a boxplot with a jittered scatter plot overlay
p_group = ggplot(df, aes(x = group, y=Percentage,fill = group))+
  facet_wrap(~cellClass, scales = "free", ncol = 3)+
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(y = Percentage, color = genotype, size = 0.5), position = position_jitter(width = 0.3, height = 0)) +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  labs(title = "Relative abundance of major cellTypes", x = "cellType", y = "Percentage")
p_group


pdf("cellClass_abundance_2groups.pdf", h =6, w = 7)
p_group
dev.off()


# Perform a two-way ANOVA
model <- aov(Ratio ~ cellClass * genotype, data = df)

# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellClass:genotype`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("PC:P1cKO-PC:CTRL","PC:P2cKO-PC:CTRL","PC:dKO-PC:CTRL",
                                               "IN:P1cKO-IN:CTRL","IN:P2cKO-IN:CTRL","IN:dKO-IN:CTRL",
                                               "CN.gaba:P1cKO-CN.gaba:CTRL","CN.gaba:P2cKO-CN.gaba:CTRL","CN.gaba:dKO-CN.gaba:CTRL",
                                               "CN:P1cKO-CN:CTRL","CN:P2cKO-CN:CTRL","CN:dKO-CN:CTRL",
                                               "GC:P1cKO-GC:CTRL","GC:P2cKO-GC:CTRL","GC:dKO-GC:CTRL",
                                               "NPC:P1cKO-NPC:CTRL","NPC:P2cKO-NPC:CTRL","NPC:dKO-NPC:CTRL"),]


filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
