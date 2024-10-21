
library(Seurat);library(tidyverse)

# reload data
# cdS <- readRDS("data/integrated_all.rds")
cdS = readRDS("/Volumes/SSD_JamesLi/Experiments1/Multiplexed/Integrated_all/data/integrated_all.rds")

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

# Examine the relative abundance of PC subtypes
meta_pc = meta %>% 
  filter(!ID1 %in% c("P1_X44675","P2_X44675","P1_X44677","P1_X49530","P2_X49531")) %>% 
  filter(startsWith(cellTypeAllM, "PC"))

tab_pc = prop.table(table(meta_pc$cellTypeAllM, meta_pc$sample), margin = 2) %>% as.data.frame()

# set column names for the dataframe
colnames(tab_pc) <- c("cellType", "sample", "Ratio")
# tab_pc$cellType = factor(tab_pc$cellType, levels = paste0("PC",1:11))

# add a new column to the dataframe to show the percentage of the cluster in the sample
tab_pc$Percentage = tab_pc$Ratio * 100

forGroup = meta_pc[!duplicated(meta_pc$sample),]
forGroup$genotype1 = forGroup$genotype
forGroup$genotype1[forGroup$genotype1 %in% c("P2cKO", "dKO")] = "P2.deficient"

forGroup$group = forGroup$genotype1
forGroup$group[forGroup$group == "P1cKO"] = "CTRL"
forGroup$group1 = forGroup$genotype
forGroup$group1[forGroup$group1 == "P1cKO"] = "CTRL"

# create a new column called genotype and add it to the dataframe
tab_pc$genotype = plyr::mapvalues(tab_pc$sample, forGroup$sample,forGroup$genotype)
tab_pc$genotype1 = plyr::mapvalues(tab_pc$sample, forGroup$sample,forGroup$genotype1)
tab_pc$group = plyr::mapvalues(tab_pc$sample, forGroup$sample,forGroup$group)
tab_pc$group1 = plyr::mapvalues(tab_pc$sample, forGroup$sample,forGroup$group1)
tab_pc$group1 = factor(tab_pc$group1, levels = c("CTRL","P2cKO","dKO"))
tab_pc$genotype = factor(tab_pc$genotype, levels = c("CTRL","P1cKO","P2cKO","dKO"))
tab_pc$genotype1 = factor(tab_pc$genotype1, levels = c("CTRL","P1cKO","P2.deficient"))

# sort PC subtype based on the changes due to the loss of Foxp2
res = aov(Ratio~group*cellType, tab_pc)
summary(res)
res1 = TukeyHSD(res)$`group:cellType` %>% as.data.frame()
res1 = res1[c("P2.deficient:PC5-CTRL:PC5","P2.deficient:PC8-CTRL:PC8",
              "P2.deficient:PC2-CTRL:PC2","P2.deficient:PC4-CTRL:PC4",
              "P2.deficient:PC9-CTRL:PC9","P2.deficient:PC10-CTRL:PC10",
              "P2.deficient:PC6-CTRL:PC6","P2.deficient:PC7-CTRL:PC7",
              "P2.deficient:PC11-CTRL:PC11","P2.deficient:PC3-CTRL:PC3",
              "P2.deficient:PC1-CTRL:PC1"),]
res1 = res1[order(res1$`p adj`),]
order = gsub(".*:","", rownames(res1))
tab_pc$cellType = factor(tab_pc$cellType, levels = order)

p_group = ggplot(tab_pc, aes(x=group, y = Percentage, fill = group))+
  facet_wrap(~cellType, scales = "free", ncol = 4)+
  cowplot::theme_cowplot()+
  geom_boxplot()+
  geom_point(aes(y = Percentage, color = genotype, size = 0.5), position = position_jitter(width = 0.3, height = 0)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())+
  xlab("")+ylab("Percentage over total PC")
p_group


p_group1 = ggplot(tab_pc, aes(x = group1, y = Percentage, fill = group1)) +
  geom_boxplot() +
  geom_point(aes(y = Percentage, color = genotype, size = 0.5), position = position_jitter(width = 0.3, height = 0)) +
  facet_wrap(~ cellType, scales = "free") +
  # Set the color palette and legend
  cowplot::theme_cowplot()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank())+
  xlab("")+ylab("Percentage over total PC")
p_group1

pdf("CellTypes_abundance_PC_2groups.pdf", w = 10, h =6)
p_group
dev.off()

pdf("CellTypes_abundance_PC_3groups.pdf", w = 10, h =6)
p_group1
dev.off()






