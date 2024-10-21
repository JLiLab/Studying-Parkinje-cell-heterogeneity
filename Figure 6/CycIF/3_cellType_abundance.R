library(Seurat);library(tidyverse)
cdS <- readRDS("data/integrated_all.rds")
# reload data
metaDD <- cdS@meta.data
# remove the posterior-most sections and non-specific cells 
#(similar results obtained without filtering)
meta_filtered = metaDD %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(!cellTypeAllM %in% c("to.filter"))

# groups cell types into major cell groups
meta_filtered$cellTypeAllM = as.character(meta_filtered$cellTypeAllM)
table(meta_filtered$cellTypeAllM)

meta_filtered$cellGroup = meta_filtered$cellTypeAllM
meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("Dentate","Fastigial","Interposed")] <- "CN"
meta_filtered$cellGroup[grep("^GCP", meta_filtered$cellGroup)]  <- "GCP"
#meta_filtered$cellGroup[meta_filtered$cellGroup %in% c("IN","GABA.Pre")] <- "IN"
meta_filtered$cellGroup[grep("Mb", meta_filtered$cellGroup)] <- "Mb"
meta_filtered$cellGroup[grep("^PC", meta_filtered$cellGroup)] <- "PC"
table(meta_filtered$cellGroup)
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


# Create a boxplot with a jittered scatter plot overlay
p2 = ggplot(df, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, scales = "free", ncol = 4) +
  geom_jitter(width = 0.25, size = 4, alpha = 0.6)+
  scale_fill_discrete(name = "cellType") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "Relative abundance of cellType", x = "cellType", y = "Percentage")
p2


p1 <- ggplot(df, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = genotype), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free", ncol =4) +
  # Set the color palette and legend
  scale_fill_discrete(name = "cellType") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "top")+
  # Add axis and title labels
  labs(title = "Relative abundance of cellType", x = "cellType", y = "Percentage")
p1

pdf("figures/All_cellTypes_abundance.pdf", h =8, w =8)
p1
p2
dev.off()


# # remove the most posterior sections, which were not included in the CTRL and Foxp1-cKO dataset
meta_pc = metaDD %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) 
#%>% filter(!is.na(PC))

tab = prop.table(table(meta_pc$cellTypeAllM, meta_pc$sample), margin = 2) %>% as.data.frame() %>% 
  filter(str_starts(Var1, "^PC"))

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
p = ggplot(tab, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = sample, size = 1), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free") +
  # Set the color palette and legend
  scale_fill_discrete(name = "cellType") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # Add axis and title labels
  labs(title = "Relative abundance of PC subtypes", x = "cellType", y = "Percentage")
p

p1 = ggplot(tab, aes(x = genotype, y = Percentage, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, scales = "free") +
  geom_jitter(width = 0.25, size = 4, alpha = 0.6)+
  scale_fill_discrete(name = "cellType") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # Add axis and title labels
  labs(title = "Relative abundance of PC subtypes", x = "cellType", y = "Percentage")
p1

p2 = ggplot(tab, aes(x = genotype1, y = Percentage, fill = genotype1)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ cellType, scales = "free") +
  geom_jitter(width = 0.25, size = 4, alpha = 0.6)+
  scale_fill_discrete(name = "cellType") +
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  # Add axis and title labels
  labs(title = "Relative abundance of PC subtypes", x = "cellType", y = "Percentage")
p2

pdf("figures/PC_subtype_abundance.pdf", h =8, w =8)
p1
p2
p
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
# diff         lwr           upr        p adj
# PC1:P2cKO-PC1:CTRL -0.02001592 -0.03994563 -8.620728e-05 4.740897e-02
# PC1:dKO-PC1:CTRL   -0.03255247 -0.05248219 -1.262276e-02 2.033281e-06

# Perform a two-way ANOVA
model <- aov(Ratio ~ cellType * genotype1, data = tab)
# View the ANOVA table
summary(model)

# Conduct a post hoc test (Tukey's HSD)
posthoc <- TukeyHSD(model)

# Filter and show specific pairwise comparisons
specific_comparisons <- posthoc$`cellType:genotype1`

# Subset specific comparisons for PC1:mouse - PC1:human, PC2:mouse - PC2:human, etc.
filtered_comparisons <- specific_comparisons[c("PC1:CTRL-PC1:P1cKO","PC1:P2.deficient-PC1:CTRL",
                                               "PC2:CTRL-PC2:P1cKO","PC2:P2.deficient-PC2:CTRL",
                                               "PC3:CTRL-PC3:P1cKO","PC3:P2.deficient-PC3:CTRL",
                                               "PC4:CTRL-PC4:P1cKO","PC4:P2.deficient-PC4:CTRL",
                                               "PC5:CTRL-PC5:P1cKO","PC5:P2.deficient-PC5:CTRL",
                                               "PC6:CTRL-PC6:P1cKO","PC6:P2.deficient-PC6:CTRL",
                                               "PC7:CTRL-PC7:P1cKO","PC7:P2.deficient-PC7:CTRL",
                                               "PC8:CTRL-PC8:P1cKO","PC8:P2.deficient-PC8:CTRL",
                                               "PC9:CTRL-PC9:P1cKO","PC9:P2.deficient-PC9:CTRL",
                                               "PC11:CTRL-PC11:P1cKO","PC11:P2.deficient-PC11:CTRL"    ),]

filtered_comparisons=as.data.frame(filtered_comparisons)
filtered_comparisons[filtered_comparisons$`p adj`<0.05,]
# diff         lwr          upr        p adj
# PC1:P2.deficient-PC1:CTRL -0.02628420 -0.04441366 -0.008154731 6.705718e-05
# PC7:P2.deficient-PC7:CTRL -0.01922005 -0.03734951 -0.001090580 2.387005e-02

