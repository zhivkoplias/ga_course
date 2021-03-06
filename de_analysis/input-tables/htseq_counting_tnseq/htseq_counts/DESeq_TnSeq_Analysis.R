source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

directory <- "/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/htseq_counts_tnseq/htseq_counts/"

sampleFiles_Serum <- grep("Serum",list.files(directory),value=TRUE)
sampleCondition_Serum <- sub("(.*Serum).*","\\1",sampleFiles_Serum)
sampleTable_Serum <- data.frame(sampleName = sampleFiles_Serum,
                          fileName = sampleFiles_Serum,
                          condition = sampleCondition_Serum)

sampleFiles_Human <- grep("HeatInactivated",list.files(directory),value=TRUE)
sampleCondition_Human <- sub("(.*HeatInactivated).*","\\1",sampleFiles_Human)
sampleTable_Human <- data.frame(sampleName = sampleFiles_Human,
                                fileName = sampleFiles_Human,
                                condition = sampleCondition_Human)

sampleFiles_BH <- grep("BH",list.files(directory),value=TRUE)
sampleCondition_BH <- sub("(.*BH).*","\\1",sampleFiles_BH)
sampleTable_BH <- data.frame(sampleName = sampleFiles_BH,
                                fileName = sampleFiles_BH,
                                condition = sampleCondition_BH)

sampleCondition <- c(sampleCondition_Serum, sampleCondition_BH)

sampleTable <- merge(sampleTable_BH, sampleTable_Serum, all=TRUE)
#sampleTable <- merge(sampleTable_SH, sampleTable_Human, all=TRUE)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#prefiltering

keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- ddsHTSeq 

#Note on factor levels

#dds$condition <- factor(dds$condition, levels = c("Serum","Human"))
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
#biocLite("apeglm")
library(apeglm)

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm")
resOrdered <- resLFC[order(resLFC$padj),]
resOrdered

# p-values ordering

resOrdered <- res[order(res$log2FoldChange),]
summary(res)

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.1, na.rm=TRUE)

#plots DE

plotMA(res, ylim=c(-9,9))

plotMA(resLFC, ylim=c(-10,10), main='MA plot', alpha=0.1)


idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# plots counts

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.ontop = TRUE
  )

# pca

# install.packages("factoextra")
library(factoextra)

setwd("/home/erik/sweden/courses/2nd_semester/Genome_Analysis/ga_course/htseq_counts_tnseq/htseq_counts/for_pca")
temp = list.files(pattern="*.txt")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i], header=TRUE, sep="\t"))

BH1 <- merge(BHI_ERR1801012_transcripts.txt, BHI_ERR1801013_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
BH2 <- merge(BH1, BHI_ERR1801014_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H1 <- merge(BH2, Human_ERR1801009_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H2 <- merge(H1, Human_ERR1801010_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
H3 <- merge(H2, Human_ERR1801011_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S1 <- merge(H3, Serum_ERR1801006_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S2 <- merge(S1, Serum_ERR1801007_transcripts.txt, by.x="aaaGene", by.y="aaaGene")
S3 <- merge(S2, Serum_ERR1801008_transcripts.txt, by.x="aaaGene", by.y="aaaGene")

dataset_for_pca <- S3
dataset_for_pca[, 2:10] <- sapply(dataset_for_pca[, 2:10], as.numeric)
dataset_for_pca <- dataset_for_pca[, 2:10]
str(dataset_for_pca)
#install.packages("ggfortify")

res.pca = prcomp(t(as.matrix(dataset_for_pca)))
project.pca <- res.pca
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
# attemot 1
project.pca$x
autoplot(project.pca)

require(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)

scatterplot3d(project.pca$x[,1:3], angle=-40, main="", color="black", pch=17, xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), zlab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"), grid=FALSE, box=FALSE)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(project.pca$x[,1:3], grid = c("xy", "xz", "yz"))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(project.pca$x[,1:3], grid = c("xy", "xz", "yz"))

# PCA with colour

library(FactoMineR)
iris.pca <- PCA(iris, quali.sup=5)
plot(iris.pca, habillage = 5, 
     col.hab = c("green", "blue", "red"), 
     title = "Dataset projected onto PC1-2 Subspace")

#Another attempt

category = c( "BHI","BHI","BHI","Human","Human","Human","Serum","Serum","Serum")
qplot(project.pca$x[,1],project.pca$x[,2], colour=category, size=9, xlab="PCA 1", ylab="PCA 2") 

#Third
df_out <- as.data.frame(project.pca$x)
df_out$group <- sapply( strsplit(as.character(colnames(dataset_for_pca)), "_"), "[[", 1 )
head(df_out)

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()
p

# Claudios - on DE data

rld <- rlogTransformation(dds, blind=TRUE)
#plot(0, 0, main = "", xlab = "", ylab = "", xlim = c(-6, 6), ylim = c(-10, 10), col = "white", asp = 1)
plotPCA(rld)

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
    geom_point(size = 5) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}

new_plot(rld)

# more results

mcols(res)$description


#save data

save_it <- as.data.frame(only_sign_dds_rlog)

write.csv(as.data.frame(save_it[1]), 
          file="exclude_two_sample_3_cond_down.csv")

new <- read.csv("exclude_two_sample_3_cond_down.csv", sep = ',')
names <- new$X
write.csv(names, 
          file="exclude_two_sample_3_cond_gene_names_down.csv", row.names = FALSE)

#kalles code

# Create the report from the DESeq2 package, 
DESeq2Report(dds = dds, intgroup = "growth_medium")

# Since expression data is 
rld <- rlog(dds)
resultsdds <- results(dds)

only_sign_dds_rlog <- subset(resultsdds, padj < 0.1)
only_sign_dds_rlog <- subset(only_sign_dds_rlog, log2FoldChange < -2)


change_filter <- change_filter <- which(
  abs(only_sign_dds_rlog$padj < 0.1) | 
    abs(only_sign_dds_rlog$padj < 0.1)
)

topVarGenes <- head(order(-rowVars(assay(rld))),40)
#topVarGenes <- topVarGenes[(topVarGenes %in% only_sign_dds_rlog$)]

mat = assay(rld)[ head(order(res$padj),40), ] # select the top 30 genes with the lowest padj
mat = mat - rowMeans(mat) # Subtract the row means from each value
df <- as.data.frame(colData(rld)[,"condition"])
colnames(mat) <- c("BH-1","BH-2","BH-3","Serum-1","Serum-2","Serum-3")
colnames(df) <- "Condition"
rownames(df) <- colnames(mat)
#install.packages("pheatmap")
library(pheatmap)
pheatmap(mat, annotation_col=df)




# Example of DESeq analysis
plotDispEsts(dds)
log_dds<-rlog(dds)

counts_table = counts( dds, normalized=TRUE )

filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
head(filtered_norm_counts)

res_ordered<-res[order(res$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_genes<-cbind(res_ordered,GeneID)
dim(res_genes)
head(res_genes)
dim(res_genes)
res_genes_merged <- merge(res_genes,filtered_norm_counts,by=unique("GeneID"))
dim(res_genes_merged)
head(res_genes_merged)
res_ordered<-res_genes_merged[order(res_genes_merged$padj),]

plot(log2(res_ordered$baseMean), res_ordered$log2FoldChange, col=ifelse(res_ordered$padj < 0.1, "red","gray67"),main="MA plot (padj<0.05)", xlab = "Mean of normalized counts", ylab = "Log fold change", xlim=c(1,20),pch=20,cex=1,ylim=c(-5,5))
abline(h=c(-1,1), col="blue")

resSig = res_ordered[res_ordered$padj < 0.1, ]
resSig = resSig[resSig$log2FoldChange > 1 | resSig$log2FoldChange < -1,]

genes<-resSig$GeneID
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=genes,pos=1,cex=0.70)

#Gene clustering
library("genefilter")
vsd <- rlog(dds, blind = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("BH","Serum")])
pheatmap(mat, annotation_col = df)

