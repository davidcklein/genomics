setwd("~/Directory_containing_counts_files")
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
#design should include all conditions affecting the experiments; in this case, batch/replicate, gene knockout, and treatment). As an example, this script will use BRD9 KO vs ASCR (no KO) for the comparison of interest.  
dds <- DESeq(dds)
resultsNames(dds) # lists the comparisons - ensure that these are correct
res <- results(dds, contrast=c("KO","BRD9","ASCR"))
res05 <- results(dds, alpha=0.05, contrast=c("KO","BRD9","ASCR"))
#Setting a new alpha value will change the significance cutoff from the default 0.1 to 0.05
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
reslfc <- lfcShrink(dds, coef=c("KO_BRD9_vs_ASCR"), type="apeglm")
#lfcShrinkage will preserve important differences while shrinking log2-fold changes to minimize the effects of lowly expressed genes. 
reslfc
summary(reslfc,alpha=0.05, contrast=c("KO","BRD9","ASCR"))
sum(reslfc$padj < 0.05, na.rm=TRUE)
head(res)
resOrdered <- reslfc[order(reslfc$padj),]
resOrdered
plotCounts(dds, gene=which.min(reslfc$padj), intgroup="KO")
#This step is a nice sanity check - it will plot the counts of the most significantly different gene (using KO as the comparison of interest)
write.csv(as.data.frame(resOrdered), 
          file="BRD9_KO_vs_noKO_DESeq2_results.csv")
#This step will write a DESeq2 output table as a comma-delimited file sorted by ascending adjusted p-value
options(ggrepel.max.overlaps = 5)
library(ggplot2)
library("EnhancedVolcano")
#THIS COMMAND WILL PROBABLY NEED TO BE TAILORED TO THE INDIVIDUAL EXPERIMENT
#For a nice volcano plot:
  #adjust the xlim values to include the highest and lowest LFC values
  #adjust the ylim values to include the lowest adj. p value
  #To call out specific genes, add the following line:
#    selectLab = c('gene1','gene2'),
EnhancedVolcano(resOrdered,
                lab = (row.names(resOrdered)),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'BRD9 KO vs no KO',
                xlim = c(-10,10),
                ylim = c(0,-log10(10e-250)),
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

# The remainder of the script is focused on generating principal component analysis (PCA) plots. The file names are specific to the individual experiment that this script was modeled on, but should be adjusted to match the order of files in the dds data matrix.

library(ggplot2)
library(ggrepel)
rld <- rlog(dds,blind=FALSE)
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)
plotPCA(rld,intgroup=c("KO"))
PCA<-plotPCA(rld,intgroup=c("KO"))
df.name<-(c("BAF180_Cu","BAF180_Cu","BAF180_NT","BAF180_NT",
            "BAF180_Zn","BAF180_Zn","BAF250_Cu","BAF250_Cu","BAF250_NT","BAF250_NT", "BAF250_Zn","BAF250_Zn","BRD9_Cu","BRD9_Cu","BRD9_NT","BRD9_NT", "BRD9_Zn","BRD9_Zn","SCR_Cu","SCR_Cu","SCR_NT","SCR_NT", "SCR_Zn","SCR_Zn","WT_Cu","WT_Cu","WT_NT","WT_NT", "WT_Zn","WT_Zn"))
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

#This last bit will just write a table including the PC1 and PC2 values. As long as your file name matrix is properly annotated, it won't do much for you, but you can confirm that files are labeled appropriately by matching the PC values to the PCA plot points.
pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData
