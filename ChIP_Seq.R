# Human annotation using CHIPseeker 

setwd("F:\\\\zhangdw\\\\Lenovo computer\\\\20190716 backup\\\\Bioinformatics\\\\Project\\\\Trisomy 21\\\\database of birth defect\\\\DNA_protein_interaction\\\\GSE55506\\\\annotation")
library("ChIPseeker")
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library("org.Hs.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("org.Mm.eg.db")
T21 <- readPeakFile("T21_peaks_chrID.narrowPeak")
D21 <- readPeakFile("D21_peaks_chrID.narrowPeak")
peaks <- list(T21=T21, D21=D21)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnnoList <- lapply(peaks, annotatePeak, tssRegion=c(-2000,2000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000, annoDb="org.Hs.eg.db")
as.GRanges(peakAnnoList[[1]])
as.data.frame(peakAnnoList[[1]])
write.table(as.data.frame(peakAnnoList[[1]]),file = "Trisomy_21_annotation.txt",col.names = TRUE,row.names = FALSE, sep = "\\t")
write.table(as.data.frame(peakAnnoList[[2]]),file = "Euploid_annotation.txt",col.names = TRUE,row.names = FALSE, sep = "\\t")


Ts1Rhr <- readPeakFile("Ts1Rhr_H3K27me3_peaks_chr.narrowPeak")
WT <- readPeakFile("WT_H3K27me3_peaks_chr.narrowPeak")
peaks <- list(Ts1Rhr=Ts1Rhr, WT=WT)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peakAnnoList <- lapply(peaks, annotatePeak, tssRegion=c(-2000,2000), TxDb=txdb, addFlankGeneInfo=TRUE, flankDistance=5000, annoDb="org.Mm.eg.db")
as.GRanges(peakAnnoList[[1]])
write.table(as.data.frame(peakAnnoList[[1]]),file = "Ts1Rhr_H3K27ac_annotation.txt",col.names = TRUE,row.names = FALSE, sep = "\t")
write.table(as.data.frame(peakAnnoList[[2]]),file = "WT_H3K27ac_annotation.txt",col.names = TRUE,row.names = FALSE, sep = "\t")


# DiffBind
setwd("F:\\zhangdw\\Lenovo computer\\20190716 backup\\Bioinformatics\\Project\\Trisomy 21\\database of birth defect\\DNA_protein_interaction\\GSE55506\\Differential_expression")
library(DiffBind)
dbObj <- dba(sampleSheet="SampleSheet.csv")
# jisuan qin he ju  zheng
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
# differential experission
# Establishing a contrast
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, , minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# summary of results
dba.show(dbObj, bContrasts=T)
# overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
# extract the result
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# save the result
# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="Differential_peaks_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="Differential_peaks_deseq2.txt", sep="\t", quote=F, col.names = NA)
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="Differential_peaks_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="Differential_peaks_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)
