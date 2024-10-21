library(Seurat);library(COSG);library(clusterProfiler);library(org.Mm.eg.db)
library(tidyverse)

# PCin <- readRDS("/Users/jamesli/Manuscripts/PC_Manuscript/data/PCin_new.rds")
PCin <- readRDS("/Users/jamesli/Manuscripts/PC_Manuscript/data_analysis/data/PCin_new.rds")
Idents(PCin)=PCin$PC0.8

mGenes = readRDS("~/Desktop/Annotation1.rds")
marker_cosg <- cosg(PCin, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=20)
markers<- marker_cosg$names %>% unclass() %>%  stack()
scores <- marker_cosg$scores %>% unclass() %>% stack()
df <- data.frame(gene = markers$values, score = scores$values, cluster = scores$ind)
df$TF = mGenes$Type[match(df$gene, mGenes$mgi_symbol)]
df$Entrez = mGenes$entrezgene[match(df$gene, mGenes$mgi_symbol)]
df$description = mGenes$description[match(df$gene, mGenes$mgi_symbol)]
openxlsx::write.xlsx(df, "tables/PC_markers_16.xlsx")

marker_cosg <- cosg(integrated, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=4)
marker_cosg <- cosg(PCin, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=3)
all_markers<- marker_cosg$names %>%
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

p <- DotPlot(object = PCin, features = all_markers)+ coord_flip()+ RotatedAxis()
df <- data.frame(Gene = p$data$features.plot, avg.exp = p$data$avg.exp.scaled, pct.exp = p$data$pct.exp, cluster = p$data$id)

p1 <- df %>% 
  filter(avg.exp > 0, pct.exp > 1) %>% 
  ggplot(aes(x=cluster, y = Gene, color = avg.exp, size = pct.exp)) + 
  geom_point() +
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab("")+
  theme(axis.ticks = element_blank()) 

pdf("../figures/PC_markers.pdf", h =6, w=5)
p1
dev.off()


df1 = df %>% filter(!duplicated(df$gene), score > 0.2)
dim(df1)

Symbol <- df

bg_ids <- unique(rownames(PCin));length(bg_ids)

ego = enrichGO(gene = df1$gene, 
                     universe      = bg_ids,
                     keyType = "SYMBOL",
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
barplot(ego, showCategory=15) 

kk <- enrichKEGG(gene         = df1$Entrez,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
barplot(kk, showCategory=15) 

mkk <- enrichMKEGG(gene = df1$Entrez,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
barplot(mkk, showCategory=15) 


mGenes <- readRDS(file="~/Desktop/Annotation1.rds")
sig_genes <- FindAllMarkers(PCin, assay = "SCT" ,logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
                            min.diff.pct = -Inf, node = NULL, verbose = TRUE,
                            only.pos = TRUE, max.cells.per.ident = Inf, random.seed = 1,
                            latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3,
                            pseudocount.use = 1, return.thresh = 0.01)
sig_genes <- mutate(sig_genes,pct.diff=pct.1-pct.2)
sig_genes$TF <- "no"
id <- intersect(mGenes$mgi_symbol,sig_genes$gene)
x <- sig_genes[sig_genes$gene %in% id,]
sig_genes$TF[sig_genes$gene %in% id] <- mGenes$Type[match(x$gene,mGenes$mgi_symbol)]
sig_genes$Entrez <- mGenes$entrezgene[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes$description <- mGenes$description[match(as.character(sig_genes$gene),mGenes$mgi_symbol)]
sig_genes <- sig_genes[, c("gene","TF","Entrez","description","p_val_adj","avg_log2FC","pct.1","pct.2","pct.diff","cluster")]
openxlsx::write.xlsx(sig_genes, "tables/PC_markers_embryonic.xlsx", overwrite = T)

res=sig_genes %>% filter(avg_log2FC>.4)
table(res$cluster)
length(unique(res$gene))

counts = GetAssayData(PCin, slot = "counts",assay = "SCT")
dim(counts)

ego = enrichGO(gene = res$gene, 
               universe      = rownames(counts),
               keyType = "SYMBOL",
               OrgDb         = org.Mm.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)
barplot(ego, showCategory=10) 

pdf("figures/enrichGO_PC.pdf", h =5, w =6)
barplot(ego, showCategory=10) 
dev.off()

kk <- enrichKEGG(gene         = res$Entrez,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
barplot(kk, showCategory=15) 

mkk <- enrichMKEGG(gene = res$Entrez,
                   organism = 'mmu',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
barplot(mkk, showCategory=15) 

write.csv(ego@result, "ego_results.csv")

ego@result$geneID[[1]]

gene_string <- ego@result$geneID[[1]]
gene_list <- strsplit(gene_string, split = "/")[[1]]
