### Load libraries

source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)

### Load data

directory <- "/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/input-tables/htseq_counting/for-r"
setwd("/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/input-tables/htseq_counting/for-r")

sampleFiles_Serum <- grep("Serum",list.files(directory),value=TRUE)
sampleCondition_Serum <- sub("(.*Serum).*","\\1",sampleFiles_Serum)
sampleTable_Serum <- data.frame(sampleName = sampleFiles_Serum,
                          fileName = sampleFiles_Serum,
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

### Prefiltering

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- ddsHTSeq 

### Choosing the reference condition

dds$condition <- relevel(dds$condition, ref = "BH")

### DE analysis

dds <- dds[, dds$condition %in% c("Serum","BH") ]
dds <- DESeq(dds)
res <- results(dds)
res

res <- results(dds, name="condition_Serum_vs_BH")
res <- results(dds, contrast=c("condition","Serum","BH"))

### Plot PCA

rld <- rlogTransformation(dds, blind=TRUE)

new_plot <- function (object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    ylim(-20, 20) +
    theme_bw() +
    geom_point(size = 10) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                         100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                             100), "% variance")) + coord_fixed()
}

new_plot(rld)


### Plot heatmap with variance of LogFoldChange

rld <- rlog(dds)
resultsdds <- results(dds)

topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat = assay(rld)[ head(order(res$padj),20), ] # select the top 20 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value 
df <- as.data.frame(colData(rld)[,"condition"])
colnames(mat) <- c("BH-1","BH-2","BH-3","Serum-1","Serum-2","Serum-3")
colnames(df) <- "Condition"
rownames(df) <- colnames(mat)

pheatmap(mat, annotation_col=df)

### p-values ordering

resOrdered <- res[order(-res$log2FoldChange),]
summary(res)

sum(res$padj < 0.001, na.rm=TRUE)

res05 <- results(dds, alpha=0.001)
summary(res05)

sum(res05$padj < 0.001, na.rm=TRUE)


### plots DE with MA

plotMA(res05, ylim=c(-9,9))

### Save data

save_it <- as.data.frame(resOrdered)
save_it <- subset(save_it, padj<0.001)
save_it_up <- subset(save_it, log2FoldChange > 2 )
save_it_down <- subset(save_it, log2FoldChange < 0.5 )

write.table(save_it_up, 
            file="/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/rna-seq/genes_up_in_Serum.csv",row.names=TRUE, sep='\t')
write.table(save_it_up[0], 
          file="/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/rna-seq/genes_up_in_Serum_gene_names.csv", row.names=TRUE, sep='\t')
          

write.table(save_it_down, 
            file="/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/rna-seq/genes_down_in_Serum.csv",row.names=TRUE, sep='\t')
write.table(save_it_down[0], 
            file="/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/de_analysis/rna-seq/genes_down_in_Serum_gene_names.csv", row.names=TRUE, sep='\t')