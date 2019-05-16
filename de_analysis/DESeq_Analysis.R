#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

directory <- "/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/htseq_counting"

sampleFiles_Serum <- grep("Serum",list.files(directory),value=TRUE)
sampleCondition_Serum <- sub("(.*Serum).*","\\1",sampleFiles_Serum)
sampleTable_Serum <- data.frame(sampleName = sampleFiles_Serum,
                          fileName = sampleFiles,
                          condition = sampleCondition_Serum)

sampleFiles_BH <- grep("BH",list.files(directory),value=TRUE)
sampleCondition_BH <- sub("(.*BH).*","\\1",sampleFiles_BH)
sampleTable_BH <- data.frame(sampleName = sampleFiles_BH,
                                fileName = sampleFiles_BH,
                                condition = sampleCondition_BH)

sampleCondition <- c(sampleCondition_Serum, sampleCondition_BH)
sampleTable <- merge(sampleTable_BH, sampleTable_Serum, all=TRUE)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#prefiltering

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- ddsHTSeq 

#Note on factor levels

#dds$condition <- factor(dds$condition, levels = c("Serum","BH"))
#dds$condition <- relevel(dds$condition, ref = "Serum")
#dds$condition <- droplevels(dds$condition)

# DE analysis
dds <- dds[, dds$condition %in% c("Serum","BH") ]
dds <- DESeq(dds)
res <- results(dds)
res

res <- results(dds, name="condition_Serum_vs_BH")
res <- results(dds, contrast=c("condition","Serum","BH"))


# Log fold change shrinkage for visualization and ranking
biocLite("apeglm")

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm")
resLFC

# p-values ordering

resOrdered <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

#plots DE

plotMA(res, ylim=c(-9,9))

plotMA(resLFC, ylim=c(-9,9))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# plots counts

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# more results

mcols(res)$description


#save data

write.csv(as.data.frame(resOrdered), 
          file="condition_Serum_vs_BH_results.csv")
