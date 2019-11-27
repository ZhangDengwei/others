library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)

up <- read.csv("up.csv", header = T)
gene <- up$Gene_name
gene <- gsub(" ","",gene)
test <- bitr(gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),OrgDb = "org.Mm.eg.db")
ego <- enrichGO(gene = test$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL")
kk <- enrichKEGG(gene = test$ENTREZID, organism = 'mmu', pvalueCutoff = 1)

pdf("kkup.pdf", 13,10)
barplot(kk,showCategory = 20) + ggtitle("barplot for kegg")
dotplot(kk, showCategory=20) + ggtitle("dotplot for kegg")
try(kk <- setReadable(kk, 'org.Mm.eg.db', 'ENTREZID'))
try(cnetplot(kk))
dev.off()

pdf("goup.pdf", 13,10)
barplot(ego,showCategory = 20) + ggtitle("barplot for go")
dotplot(ego, showCategory=20) + ggtitle("dotplot for go")
try(ego <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID'))
try(cnetplot(ego))
dev.off()

write.table(kk, file = "kkup.txt", sep = "\t")
write.table(ego, file = "goup.txt", sep = "\t")

################################################################################################

down <- read.csv("down.csv", header = T)
gene2 <- down$Gene_name
gene2 <- gsub(" ","",gene2)
test2 <- bitr(gene2, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),OrgDb = "org.Mm.eg.db")
ego_2 <- enrichGO(gene = test2$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL")
kk2 <- enrichKEGG(gene = test2$ENTREZID, organism = 'mmu', pvalueCutoff = 1)

pdf("kkdown.pdf", 15,10)
barplot(kk2,showCategory = 20) + ggtitle("barplot for kegg")
dotplot(kk2, showCategory=20) + ggtitle("dotplot for kegg")
try(kk2 <- setReadable(kk2, 'org.Mm.eg.db', 'ENTREZID'))
try(cnetplot(kk2))
dev.off()

pdf("godown.pdf", 15,10)
barplot(ego_2,showCategory = 20) + ggtitle("barplot for go")
dotplot(ego_2, showCategory=20) + ggtitle("dotplot for go")
try(oragnx <- setReadable(ego_2, 'org.Mm.eg.db', 'ENTREZID'))
try(cnetplot(ego_2))
dev.off()

write.table(kk2, file = "kkdown.txt", sep = "\t")
write.table(ego_2, file = "godown.txt", sep = "\t")