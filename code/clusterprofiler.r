###### Code if we have gene symbols as an input ####

library(clusterProfiler)
library(ggplot2)
de <- read.csv("gene_list", header = TRUE)                       #### Give Entrez id only and not symbol. Put header as well
ego <- enrichGO(de$Entrezid, OrgDb = "org.Hs.eg.db", ont="BP", readable = TRUE, minGSSize = 10, maxGSSize = 500, keyType="SYMBOL")
ego2 <- simplify(ego, cutoff=0.8, by="p.adjust", select_fun=min, measure = "Wang")
dotplot(ego2, x = "GeneRatio", color = "p.adjust", showCategory = 20, font.size = 10, label_format = 50)
ggsave("Image.png", units="in", width=7, height=6, dpi=600)
write.csv(ego2, file = "file")

