library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)

up <- read.csv("up.csv", header = T)
gene <- as.character(up$Symbol)
test <- bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = test$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL")
kk <- enrichKEGG(gene = test$ENTREZID, organism = 'hsa', pvalueCutoff = 1)

pdf("kkup.pdf", 13,10)
barplot(kk,showCategory = 20) + ggtitle("barplot for kegg")
dotplot(kk, showCategory=20) + ggtitle("dotplot for kegg")
try(kk <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID'))
try(cnetplot(kk))
dev.off()

pdf("goup.pdf", 13,10)
barplot(ego,showCategory = 20) + ggtitle("barplot for go")
dotplot(ego, showCategory=20) + ggtitle("dotplot for go")
try(ego <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID'))
try(cnetplot(ego))
dev.off()

write.table(kk, file = "kkup.txt", sep = "\t")
write.table(ego, file = "goup.txt", sep = "\t")

################################################################################################

down <- read.csv("down.csv", header = T)
gene2 <- down$Symbol
test2 <- bitr(gene2, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
ego_2 <- enrichGO(gene = test2$SYMBOL,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")
kk2 <- enrichKEGG(gene = test2$ENTREZID,organism = 'hsa',pvalueCutoff = 1)

pdf("kkdown.pdf", 13,10)
barplot(kk2,showCategory = 20) + ggtitle("barplot for kegg")
dotplot(kk2, showCategory=20) + ggtitle("dotplot for kegg")
try(kk2 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID'))
try(cnetplot(kk2))
dev.off()

pdf("godown.pdf", 13,10)
barplot(ego_2,showCategory = 20) + ggtitle("barplot for go")
dotplot(ego_2, showCategory=20) + ggtitle("dotplot for go")
try(ego_2 <- setReadable(ego_2, 'org.Hs.eg.db', 'ENTREZID'))
try(cnetplot(ego_2))
dev.off()

write.table(kk2, file = "kkdown.txt", sep = "\t")
write.table(ego_2, file = "godown.txt", sep = "\t")