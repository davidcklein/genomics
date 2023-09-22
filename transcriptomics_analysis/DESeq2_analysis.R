setwd("~/Directory_containing_counts_file(s)")
library(DESeq2)
library(apeglm)
library(svglite)
cts <- read.table("counts.txt",header=TRUE,row.names=1)
head(cts)
coldata<-read.table("coldata.txt",header=TRUE)
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~batch+KO+treatment)
dds <- DESeq(dds)
resultsNames(dds) # lists the comparisons
BRD9_KO_res <- results(dds, contrast=c("KO","BRD9","ASCR"))
BRD9_KO_res05 <- results(dds, alpha=0.05, contrast=c("KO","BRD9","ASCR"))
summary(BRD9_KO_res05)
sum(BRD9_KO_res05$padj < 0.05, na.rm=TRUE)
BRD9_KO_reslfc <- lfcShrink(dds, coef=c("KO_BRD9_vs_ASCR"), type="apeglm")
BRD9_KO_reslfc
summary(BRD9_KO_reslfc,alpha=0.05, contrast=c("KO","BRD9","ASCR"))
sum(BRD9_KO_reslfc$padj < 0.05, na.rm=TRUE)
head(BRD9_KO_res)
BRD9_KO_resOrdered <- BRD9_KO_res[order(BRD9_KO_res$padj),]
BRD9_KO_resOrdered
plotCounts(dds, gene=which.min(BRD9_KO_reslfc$padj), intgroup="KO")
write.csv(as.data.frame(BRD9_KO_resOrdered), 
          file="BRD9_KO_vs_noKO_DESeq.csv")
options(ggrepel.max.overlaps = 5)
library(ggplot2)
library("EnhancedVolcano")
EnhancedVolcano(BRD9_KO_resOrdered,
                lab = (row.names(BRD9_KO_resOrdered)),
                x = 'log2FoldChange',
                y = 'padj',
                
                title = 'BRD9 KO vs no KO',
                xlim = c(-4,4),
                ylim = c(0,-log10(10e-50)),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.5,
                labSize = 4.0,
                labCol = 'black',
                col=c('grey', 'lightblue', 'royalblue', 'red'),
                colAlpha = 0.45, 
                legendLabels=c('NS','Log2FC','padj',
                               'padj and Log2FC'),
                legendPosition = 'none',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = FALSE,
                endsConnectors = 'last',
                lengthConnectors=unit(0.02,'npc'),
                widthConnectors = 0.4)
ggsave("BRD9_KO_vs_noKO_volcanoplot.svg", dpi=300, limitsize = FALSE)


library(ggplot2)
library(ggrepel)
rld <- rlog(dds,blind=FALSE)
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)
plotPCA(rld,intgroup=c("KO"))
PCA<-plotPCA(rld,intgroup=c("KO"))
df.name<-(c("BAF180_Cu","BAF180_Cu","BAF180_NT","BAF180_NT",
            "BAF180_Zn","BAF180_Zn","BAF250_Cu","BAF250_Cu","BAF250_NT","BAF250_NT",
            "BAF250_Zn","BAF250_Zn","BRD9_Cu","BRD9_Cu","BRD9_NT","BRD9_NT",
            "BRD9_Zn","BRD9_Zn","SCR_Cu","SCR_Cu","SCR_NT","SCR_NT",
            "SCR_Zn","SCR_Zn","WT_Cu","WT_Cu","WT_NT","WT_NT",
            "WT_Zn","WT_Zn"))
PCA+geom_text_repel(aes(label=df.name),
                    fontface=c("bold"),size=0,force=10,
                    min.segment.length=5)+
  theme_classic(base_size=12)+coord_equal(ratio=2)+geom_point(size=3)
ggsave("PCA_coloredbyKO_noBatch.tiff",dpi=300,limitsize=FALSE,compression="lzw",width=10,height=10)

PCA2<-plotPCA(rld,intgroup=c("KO","treatment"))
PCA2+geom_text_repel(aes(label=df.name),
                     fontface=c("bold"),size=4,force=8,
                     min.segment.length=1)+
  theme_classic(base_size=16)+coord_equal(ratio=1.5)+geom_point(size=4)
ggsave("PCA_coloredbyKO+treatment.tiff",dpi=300,limitsize=FALSE,compression="lzw",width=10,height=10)

PCA3<-plotPCA(rld,intgroup=c("treatment"))
PCA3+geom_text_repel(aes(label=df.name),
                     fontface=c("bold"),size=4,force=20,
                     min.segment.length=1)+
  theme_classic(base_size=16)+coord_equal(ratio=1.25)+geom_point(size=4)
ggsave("PCA_coloredbytreatment.tiff",dpi=300,limitsize=FALSE,compression="lzw",width=10,height=10)

assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch, rld$KO)
PCA4<-plotPCA(rld,intgroup=c("treatment"))
PCA4+geom_text_repel(aes(label=df.name),
                     fontface=c("bold"),size=4,force=8,
                     min.segment.length=1)+
  theme_classic(base_size=16)+coord_equal(ratio=1)+geom_point(size=4)
ggsave("PCA_coloredbytreatment_limma_correct_KO+rep.tiff",dpi=300,limitsize=FALSE,compression="lzw",width=10,height=10)


pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData
